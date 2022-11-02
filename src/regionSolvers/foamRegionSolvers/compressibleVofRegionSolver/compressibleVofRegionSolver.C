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

#include "compressibleVofRegionSolver.H"
#include "noPhaseChange.H"
#include "hydrostaticInitialisation.H"
#include "addToRunTimeSelectionTable.H"

#define COMPRESSIBLE

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionSolvers
{
    defineTypeNameAndDebug(compressibleVof, 0);
    addToRunTimeSelectionTable(regionSolver, compressibleVof, dictionary);
}
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionSolvers::compressibleVof::compressibleVof
(
    dynamicFvMesh& mesh
)
:
    fluid(mesh),
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

    phaseChangePtr(twoPhaseChangeModel::New(mixture)),
    phaseChange(phaseChangePtr()),

    alpha1(mixture.alpha1()),
    alpha2(mixture.alpha2()),

    rho1(mixture.rho1()),
    rho2(mixture.rho2()),

    rho(mixture.rho()),

    p(mixture.p()),

    T(mixture.T()),

    psi1(mixture.thermo1().psi()),
    psi2(mixture.thermo2().psi()),

    pMin("pMin", dimPressure, mixture),

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

    pressureReference
    (
        p,
        p_rgh,
        pimple.dict()
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

    dgdt(alpha1*fvc::div(phi)),

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


    turbulence
    (
        rho,
        U,
        phi,
        rhoPhi,
        alphaPhi10,
        mixture
    ),

    K("K", 0.5*magSqr(U)),

    contErr
    (
        IOobject
        (
            "contError",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimDensity/dimTime, Zero)
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

        rAU = new volScalarField
        (
            IOobject
            (
                "rAU",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimTime/dimDensity, 1)
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionSolvers::compressibleVof::~compressibleVof()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::regionSolvers::compressibleVof::moveMesh(const bool finalIter)
{
    if (correctPhi)
    {
        divU = new volScalarField
        (
            "divU0",
            fvc::div(fvc::absolute(phi, U))
        );
    }

    fvModels.preUpdateMesh();

    bool changing = fluid::moveMesh(finalIter);

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


void Foam::regionSolvers::compressibleVof::solve()
{
    bool LTS = false;
    fvMesh& mesh = mesh_;
    const Time& runTime = runTime_;

    // --- Pressure-velocity PIMPLE corrector loop
    while (pimple.loop())
    {
        fvModels.correct();

        #include "alphaControls.H"
        #include "compressibleAlphaEqnSubCycle.H"

        turbulence.correctPhasePhi();

        solveU();
        solveT();

        while (pimple.correct())
        {
            solvep();
        }

        if (pimple.turbCorr())
        {
            turbulence.correct();
        }
    }
}


void Foam::regionSolvers::compressibleVof::correctPhiField()
{
    fvMesh& mesh = mesh_;
    const Time& runTime = runTime_;

    // Calculate absolute flux
    // from the mapped surface velocity
    phi = mesh_.Sf() & Uf();

    correctUphiBCs(U, phi, true);

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

    #include "continuityErrs.H"

    // Make the flux relative to the mesh motion
    fvc::makeRelative(phi, U);
}


void Foam::regionSolvers::compressibleVof::solveU()
{
    MRF.correctBoundaryVelocity(U);

    tUEqn =
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      - fvm::Sp(contErr, U)
      + MRF.DDt(rho, U)
      + turbulence.divDevTau(U)
     ==
        fvModels.source(rho, U)
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
        K = 0.5*magSqr(U);
    }
}


void Foam::regionSolvers::compressibleVof::solveT()
{
    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T) + fvm::div(rhoPhi, T) - fvm::Sp(contErr, T)
      - fvm::laplacian(turbulence.alphaEff(), T)
      + (
            fvc::div(fvc::absolute(phi, U), p)()() // - contErr/rho*p
          + (fvc::ddt(rho, K) + fvc::div(rhoPhi, K))()()
          - (U()&(fvModels.source(rho, U)&U)()) - contErr*K
        )
       *(
           alpha1()/mixture.thermo1().Cv()()
         + alpha2()/mixture.thermo2().Cv()()
        )
     ==
        fvModels.source(rho, T)
    );

    TEqn.relax();

    fvConstraints.constrain(TEqn);

    TEqn.solve();

    fvConstraints.constrain(T);

    mixture.correctThermo();
    mixture.correct();
}


void Foam::regionSolvers::compressibleVof::solvep()
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

    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

    // Cache the phase change pressure source
    fvScalarMatrix Sp_rgh(phaseChange.Sp_rgh(rho, gh, p_rgh));

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phiHbyA, U);

    tmp<fvScalarMatrix> p_rghEqnComp1;
    tmp<fvScalarMatrix> p_rghEqnComp2;
    const surfaceScalarField& alphaPhi1 = talphaPhi1();
    surfaceScalarField alphaPhi2("alphaPhi2", phi - alphaPhi1);

    if (pimple.transonic())
    {
        #include "rhofs.H"

        surfaceScalarField phid1("phid1", fvc::interpolate(psi1)*phi);
        surfaceScalarField phid2("phid2", fvc::interpolate(psi2)*phi);

        p_rghEqnComp1 =
            (
                (fvc::ddt(alpha1, rho1) + fvc::div(alphaPhi1*rho1f))/rho1
              - fvc::ddt(alpha1) - fvc::div(alphaPhi1)
              + (alpha1/rho1)
               *correction
                (
                    psi1*fvm::ddt(p_rgh)
                  + fvm::div(phid1, p_rgh) - fvm::Sp(fvc::div(phid1), p_rgh)
                )
            );

        p_rghEqnComp2 =
            (
               (fvc::ddt(alpha2, rho2) + fvc::div(alphaPhi2*rho2f))/rho2
             - fvc::ddt(alpha2) - fvc::div(alphaPhi2)
             + (alpha2/rho2)
              *correction
               (
                   psi2*fvm::ddt(p_rgh)
                 + fvm::div(phid2, p_rgh) - fvm::Sp(fvc::div(phid2), p_rgh)
               )
           );
    }
    else
    {
        #include "rhofs.H"

        p_rghEqnComp1 =
            (
                (fvc::ddt(alpha1, rho1) + fvc::div(alphaPhi1*rho1f))/rho1
              - fvc::ddt(alpha1) - fvc::div(alphaPhi1)
              + (alpha1*psi1/rho1)*correction(fvm::ddt(p_rgh))
            );

        p_rghEqnComp2 =
            (
               (fvc::ddt(alpha2, rho2) + fvc::div(alphaPhi2*rho2f))/rho2
             - fvc::ddt(alpha2) - fvc::div(alphaPhi2)
             + (alpha2*psi2/rho2)*correction(fvm::ddt(p_rgh))
            );
    }

    if (mesh.moving())
    {
        p_rghEqnComp1.ref() += fvc::div(mesh.phi())*alpha1;
        p_rghEqnComp2.ref() += fvc::div(mesh.phi())*alpha2;
    }

    p_rghEqnComp1.ref() *= pos(alpha1);
    p_rghEqnComp2.ref() *= pos(alpha2);

    p_rghEqnComp1.ref() -=
        (fvModels.source(alpha1, mixture.thermo1().rho())&rho1)/rho1;
    p_rghEqnComp2.ref() -=
        (fvModels.source(alpha2, mixture.thermo2().rho())&rho2)/rho2;

    if (pimple.transonic())
    {
        p_rghEqnComp1.ref().relax();
        p_rghEqnComp2.ref().relax();
    }

    // Cache p_rgh prior to solve for density update
    volScalarField p_rgh_0(p_rgh);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqnIncomp
        (
            fvc::div(phiHbyA) - fvm::laplacian(rAUf, p_rgh)
         == Sp_rgh
        );

        ::Foam::solve
        (
            p_rghEqnComp1() + p_rghEqnComp2() + p_rghEqnIncomp
        );

        if (pimple.finalNonOrthogonalIter())
        {
            p = max(p_rgh + (alpha1*rho1 + alpha2*rho2)*gh, pMin);
            p_rgh = p - (alpha1*rho1 + alpha2*rho2)*gh;

            dgdt =
            (
                alpha1*(p_rghEqnComp2 & p_rgh)
              - alpha2*(p_rghEqnComp1 & p_rgh)
            );

            phi = phiHbyA + p_rghEqnIncomp.flux();

            U = HbyA
              + rAU()*fvc::reconstruct((phig + p_rghEqnIncomp.flux())/rAUf);
            U.correctBoundaryConditions();
            fvConstraints.constrain(U);
        }
    }

    // Correct Uf if the mesh is moving
    fvc::correctUf(Uf, U, fvc::absolute(phi, U));

    // Update densities from change in p_rgh
    mixture.thermo1().correctRho(psi1*(p_rgh - p_rgh_0));
    mixture.thermo2().correctRho(psi2*(p_rgh - p_rgh_0));
    mixture.correct();

    // Correct p_rgh for consistency with p and the updated densities
    p_rgh = p - rho*gh;
    p_rgh.correctBoundaryConditions();

    K = 0.5*magSqr(U);
}


Foam::scalar Foam::regionSolvers::compressibleVof::CoNum() const
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


Foam::scalar Foam::regionSolvers::compressibleVof::maxCo() const
{
    return
        max
        (
            runTime_.controlDict().lookup<scalar>("maxCo"),
            runTime_.controlDict().lookup<scalar>("maxAlphaCo")
        );
}

// ************************************************************************* //
