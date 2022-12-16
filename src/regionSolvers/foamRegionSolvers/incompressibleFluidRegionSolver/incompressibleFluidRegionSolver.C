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

#include "incompressibleFluidRegionSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionSolvers
{
    defineTypeNameAndDebug(incompressibleFluid, 0);
    addToRunTimeSelectionTable(regionSolver, incompressibleFluid, dictionary);
}
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionSolvers::incompressibleFluid::incompressibleFluid
(
    dynamicFvMesh& mesh,
    const regionSolverList& regions
)
:
    fluid(mesh, regions),
    pimple_(mesh_),
    g_
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
        linearInterpolate(U) & mesh_.Sf()
    ),

    p
    (
        IOobject
        (
            "p",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    pressureReference_
    (
        p,
        pimple_.dict()
    ),

    laminarTransport_(U, phi),
    turbulence_
    (
        incompressible::momentumTransportModel::New
        (
            U,
            phi,
            laminarTransport_
        )
    ),

    cumulativeContErr(0.0),

    MRF(mesh_),
    fvModels(fvModels::New(mesh_)),
    fvConstraints(fvConstraints::New(mesh_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionSolvers::incompressibleFluid::~incompressibleFluid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::regionSolvers::incompressibleFluid::moveMesh(const IterType iter)
{
    const fvMesh& mesh = mesh_;
    const Time& runTime = runTime_;
    fvModels.preUpdateMesh();

    bool changing = fluid::moveMesh(iter);

    if (mesh_.changing())
    {
        MRF.update();

        if (mesh.moving())
        {
            // Calculate absolute flux
            // from the mapped surface velocity
            phi = mesh.Sf() & Uf();

            correctUphiBCs(U, phi, true);

            CorrectPhi
            (
                phi,
                U,
                p,
                dimensionedScalar("rAUf", dimTime, 1),
                geometricZeroField(),
                pressureReference_,
                pimple_
            );

            #include "continuityErrs.H"

            // Make the flux relative to the mesh motion
            fvc::makeRelative(phi, U);
        }

        // if (checkMeshCourantNo)
        {
            const fvMesh& mesh = mesh_;
            const Time& runTime = runTime_;
            #include "meshCourantNo.H"
        }
    }
    return changing;
}


void Foam::regionSolvers::incompressibleFluid::solve()
{
    // --- Pressure-velocity PIMPLE corrector loop
    while (pimple_.loop())
    {
        fvModels.correct();

        solveU();

        // --- Pressure corrector loop
        while (pimple_.correct())
        {
            solvep();
        }

        if (pimple_.turbCorr())
        {
            laminarTransport_.correct();
            turbulence_->correct();
        }
    }
}


void Foam::regionSolvers::incompressibleFluid::solveU()
{
    MRF.correctBoundaryVelocity(U);

    tUEqn =
    (
        fvm::ddt(U) + fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence_->divDevSigma(U)
     ==
      - g_
      + fvModels.source(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvConstraints.constrain(UEqn);

    if (pimple_.momentumPredictor())
    {
        ::Foam::solve(UEqn == -fvc::grad(p));

        fvConstraints.constrain(U);
    }
}


void Foam::regionSolvers::incompressibleFluid::solvep()
{
    const fvMesh& mesh = mesh_;
    const Time& runTime = runTime_;

    volScalarField rAU(1.0/tUEqn().A());
    volVectorField HbyA(constrainHbyA(rAU*tUEqn().H(), U, p));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + MRF.zeroFilter(fvc::interpolate(rAU)*fvc::ddtCorr(U, phi, Uf))
    );

    MRF.makeRelative(phiHbyA);

    if (p.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p);
        fvc::makeAbsolute(phiHbyA, U);
    }

    tmp<volScalarField> rAtU(rAU);

    if (pimple_.consistent())
    {
        rAtU = 1.0/max(1.0/rAU - tUEqn().H1(), 0.1/rAU);
        phiHbyA +=
            fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
        HbyA -= (rAU - rAtU())*fvc::grad(p);
    }

    if (pimple_.nCorrPiso() <= 1)
    {
        tUEqn.clear();
    }

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, U, phiHbyA, rAtU(), MRF);

    // Non-orthogonal pressure corrector loop
    while (pimple_.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
        );

        pEqn.setReference
        (
            pressureReference_.refCell(),
            pressureReference_.refValue()
        );

        pEqn.solve();

        if (pimple_.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }
    }

    #include "continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    U = HbyA - rAtU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvConstraints.constrain(U);

    // Correct Uf if the mesh is moving
    fvc::correctUf(Uf, U, phi);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);
}

Foam::scalar Foam::regionSolvers::incompressibleFluid::CoNum() const
{
    const fvMesh& mesh = mesh_;
    const Time& runTime = runTime_;
    #include "CourantNo.H"
    return CoNum;
}


Foam::scalar Foam::regionSolvers::incompressibleFluid::maxCo() const
{
    return 1.0;
}

// ************************************************************************* //
