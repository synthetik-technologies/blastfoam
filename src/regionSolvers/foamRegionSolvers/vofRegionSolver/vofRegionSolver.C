/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "vofRegionSolver.H"
#include "noPhaseChange.H"
#include "hydrostaticInitialisation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionSolvers
{
    defineTypeNameAndDebug(vof, 0);
    addToRunTimeSelectionTable(regionSolver, vof, dictionary);
}
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionSolvers::vof::vof
(
    dynamicFvMesh& mesh,
    const regionSolverList& regions
)
:
    fluid(mesh, regions),
    pimple(mesh_),

    p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    U
    (
        IOobject
        (
            "U",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    phi
    (
        IOobject
        (
            "phi",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::flux(U)
    ),

    mixture(U, phi),

    // phaseChangePtr(twoPhaseChangeModel::New(mixture_)),
    // phaseChange(phaseChangePtr()),

    alpha1(mixture.alpha1()),
    alpha2(mixture.alpha2()),
    rho1(mixture.rho1()),
    rho2(mixture.rho2()),

    rho
    (
        IOobject
        (
            "rho",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alpha1*rho1 + alpha2*rho2
    ),

    alphaRestart(false),
    alphaPhi10
    (
        IOobject
        (
            IOobject::groupName("alphaPhi0", alpha1.group()),
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(alpha1)*phi
    ),
    rhoPhi
    (
        IOobject
        (
            "rhoPhi",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(rho)*phi
    ),

    turbulence
    (
        incompressible::momentumTransportModel::New
        (
            U,
            phi,
            mixture
        )
    ),

    g
    (
        IOobject
        (
            "g",
            runTime_.constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    hRef
    (
        IOobject
        (
            "hRef",
            runTime_.constant(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar(dimLength, Zero)
    ),
    ghRef(mag(g)*hRef),
    gh("gh", (g & mesh_.C()) + ghRef),
    ghf("ghf", (g & mesh_.Cf()) + ghRef),
    p
    (
        IOobject
        (
            "p",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*gh
    ),

    pressureReference
    (
        p,
        p_rgh,
        pimple.dict()
    ),

    correctPhi(pimple.dict().lookupOrDefault("correctPhi", true)),
    cumulativeContErr(0.0),

    MRF(mesh_),
    fvModels(fvModels::New(mesh_)),
    fvConstraints(fvConstraints::New(mesh_))
{
    if (alphaPhi10.headerOk())
    {
        Info<< "Restarting alpha" << endl;
        alphaRestart = true;
    }

    rho.oldTime();
    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pressureReference.refValue()
          - getRefCellValue(p, pressureReference.refCell())
        );
        p_rgh = p - rho*gh;
    }

    mesh_.setFluxRequired(p_rgh.name());
    mesh_.setFluxRequired(alpha1.name());

    const Time& runTime = runTime_;
    if (correctPhi)
    {
        Info<< "Constructing face velocity Uf\n" << endl;

        Uf = new surfaceVectorField
        (
            IOobject
            (
                "Uf",
                runTime_.timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(U)
        );
    }

    if
    (
        !runTime_.restart()
     // || isType<twoPhaseChangeModels::noPhaseChange>(phaseChange)
    )
    {
        if (correctPhi)
        {
            rAU = new volScalarField
            (
                IOobject
                (
                    "rAU",
                    runTime_.timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimTime/dimDensity, 1)
            );

            correctUphiBCs(U, phi, true);

            CorrectPhi
            (
                phi,
                U,
                p_rgh,
                surfaceScalarField("rAUf", fvc::interpolate(rAU())),
                geometricZeroField(),
                pressureReference,
                pimple
            );
        }
        else
        {
            correctUphiBCs(U, phi, true);

            CorrectPhi
            (
                phi,
                U,
                p_rgh,
                dimensionedScalar(dimTime/rho.dimensions(), 1),
                geometricZeroField(),
                pressureReference,
                pimple
            );
        }
    }

    #include "continuityErrs.H"
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionSolvers::vof::~vof()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::regionSolvers::vof::moveMesh(const IterType iter)
{
    if
    (
        correctPhi
     // && !isType<twoPhaseChangeModels::noPhaseChange>(phaseChange)
    )
    {
        divU = new volScalarField
        (
            "divU0",
            fvc::div(fvc::absolute(phi, U))
        );
    }

    fvModels.preUpdateMesh();

    bool changing = fluid::moveMesh(iter);

    if (mesh_.changing())
    {
        if (mesh_.topoChanging())
        {
            talphaPhi1Corr0.clear();
        }

        gh = (g & mesh_.C()) + ghRef;
        ghf = (g & mesh_.Cf()) + ghRef;

        MRF.update();

        if (correctPhi)
        {
            correctPhiField();
        }

        mixture.correct();

        // if (checkMeshCourantNo)
        {
            const fvMesh& mesh = mesh_;
            const Time& runTime = runTime_;
            #include "meshCourantNo.H"
        }
    }
    divU.clear();
    return changing;
}


void Foam::regionSolvers::vof::solve()
{
    bool LTS = false;
    fvMesh& mesh = mesh_;
    const Time& runTime = runTime_;

    // --- Pressure-velocity PIMPLE corrector loop
    while (pimple.loop())
    {
        fvModels.correct();

        #include "alphaControls.H"
        #include "alphaEqnSubCycle.H"

        mixture.correct();

        solveU();

        while (pimple.correct())
        {
            solvep();
        }

        if (pimple.turbCorr())
        {
            turbulence->correct();
        }
    }
}


void Foam::regionSolvers::vof::correctPhiField()
{
    fvMesh& mesh = mesh_;
    const Time& runTime = runTime_;

    // Calculate absolute flux
    // from the mapped surface velocity
    phi = mesh_.Sf() & Uf();

    correctUphiBCs(U, phi, true);

    if (divU.valid())
    {
        CorrectPhi
        (
            phi,
            U,
            p_rgh,
            surfaceScalarField("rAUf", fvc::interpolate(rAU())),
            divU(),
            pressureReference,
            pimple
        );
    }
    else
    {
        CorrectPhi
        (
            phi,
            U,
            p_rgh,
            surfaceScalarField("rAUf", fvc::interpolate(rAU())),
            geometricZeroField(),
            pressureReference,
            pimple
        );
    }

    #include "continuityErrs.H"

    // Make the flux relative to the mesh motion
    fvc::makeRelative(phi, U);
}


void Foam::regionSolvers::vof::solveU()
{
    MRF.correctBoundaryVelocity(U);

    tUEqn =
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevTau(rho, U)
     ==
        fvModels.source(rho, U)
      // + phaseChange.SU(rho, rhoPhi, U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvConstraints.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        ::Foam::solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh_.magSf()
            )
        );

        fvConstraints.constrain(U);
    }
}



void Foam::regionSolvers::vof::solvep()
{
    fvMesh& mesh = mesh_;
    const Time& runTime = runTime_;

    if (rAU.valid())
    {
        rAU.ref() = 1.0/tUEqn().A();
    }
    else
    {
        rAU = 1.0/tUEqn().A();
    }

    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));
    volVectorField HbyA(constrainHbyA(rAU()*tUEqn().H(), U, p_rgh));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + MRF.zeroFilter(fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi, Uf))
    );
    MRF.makeRelative(phiHbyA);

    if (p_rgh.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p_rgh);
        fvc::makeAbsolute(phiHbyA, U);
    }

    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh_.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

    // Cache the phase change pressure source
    // fvScalarMatrix Sp_rgh(phaseChange.Sp_rgh(rho, gh, p_rgh));

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvc::div(phiHbyA) - fvm::laplacian(rAUf, p_rgh)
         // == Sp_rgh
        );

        p_rghEqn.setReference
        (
            pressureReference.refCell(),
            getRefCellValue(p_rgh, pressureReference.refCell())
        );

        p_rghEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + p_rghEqn.flux();

            p_rgh.relax();

            U = HbyA + rAU()*fvc::reconstruct((phig + p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvConstraints.constrain(U);
        }
    }

    #include "continuityErrs.H"

    // Correct Uf if the mesh is moving
    fvc::correctUf(Uf, U, phi);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pressureReference.refValue()
          - getRefCellValue(p, pressureReference.refCell())
        );
        p_rgh = p - rho*gh;
    }
    p.correctBoundaryConditions();

    if (!correctPhi)
    {
        rAU.clear();
    }
}

Foam::scalar Foam::regionSolvers::vof::CoNum() const
{

    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi))().primitiveField()
    );

    scalar CoNum = 0.5*gMax(sumPhi/mesh_.V().field())*runTime_.deltaTValue();

    scalar meanCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh_.V().field()))*runTime_.deltaTValue();

    Info<< mesh_.name() << ": Courant Number mean: " << meanCoNum
        << " max: " << CoNum << endl;


    // scalar alphaCoNum = 0.0;
    // scalar meanAlphaCoNum = 0.0;

    // if (mesh.nInternalFaces())
    // {
    //     scalarField sumPhi
    //     (
    //         mixture.nearInterface()().primitiveField()
    //        *fvc::surfaceSum(mag(phi))().primitiveField()
    //     );

    //     alphaCoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    //     meanAlphaCoNum =
    //         0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
    // }

    // Info<< "Interface Courant Number mean: " << meanAlphaCoNum
    //     << " max: " << alphaCoNum << endl;

    return CoNum;//max(CoNum, alphaCoNum);
}


Foam::scalar Foam::regionSolvers::vof::maxCo() const
{
    return
        min
        (
            runTime_.controlDict().lookupOrDefault
            (
                mesh_.name() + "MaxCo",
                runTime_.controlDict().lookup<scalar>("maxCo")
            ),
            runTime_.controlDict().lookupOrDefault
            (
                mesh_.name() + "MaxAlphaCo",
                runTime_.controlDict().lookup<scalar>("maxAlphaCo")
            )
        );
}

Foam::scalar Foam::regionSolvers::vof::newDeltaT() const
{
    scalar maxCo =
        runTime_.controlDict().lookupOrDefault
        (
            mesh_.name() + "MaxCo",
            runTime_.controlDict().lookup<scalar>("maxCo")
        );
    scalar maxAlphaCo =
        runTime_.controlDict().lookupOrDefault
        (
            mesh_.name() + "MaxAlphaCo",
            runTime_.controlDict().lookup<scalar>("maxAlphaCo")
        );

    // Courant number
    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi))().primitiveField()
    );

    scalar CoNum = 0.5*gMax(sumPhi/mesh_.V().field())*runTime_.deltaTValue();

    scalar meanCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh_.V().field()))*runTime_.deltaTValue();

    Info<< mesh_.name() << ": Courant Number mean: " << meanCoNum
        << " max: " << CoNum << endl;


    scalar alphaCoNum = 0.0;
    scalar meanAlphaCoNum = 0.0;

    if (mesh_.nInternalFaces())
    {
        scalarField sumPhi
        (
            mixture.nearInterface()().primitiveField()
           *fvc::surfaceSum(mag(phi))().primitiveField()
        );

        alphaCoNum = 0.5*gMax(sumPhi/mesh_.V().field())*runTime_.deltaTValue();

        meanAlphaCoNum =
            0.5*(gSum(sumPhi)/gSum(mesh_.V().field()))*runTime_.deltaTValue();
    }

    Info<< "Interface Courant Number mean: " << meanAlphaCoNum
        << " max: " << alphaCoNum << endl;

    scalar maxDeltaTFact =
        min
        (
            maxCo/(CoNum + small),
            maxAlphaCo/(alphaCoNum + small)
        );
    scalar deltaTFact =
        min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

    return deltaTFact*runTime_.deltaTValue();
}

// ************************************************************************* //
