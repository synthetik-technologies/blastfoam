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

#include "compressibleFluidRegionSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionSolvers
{
    defineTypeNameAndDebug(compressibleFluid, 0);
    addToRunTimeSelectionTable(regionSolver, compressibleFluid, dictionary);
}
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionSolvers::compressibleFluid::compressibleFluid
(
    dynamicFvMesh& mesh
)
:
    fluid(mesh),
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
    thermo_(fluidThermo::New(mesh_)),
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
        thermo_->rho()
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
        linearInterpolate(rho*U) & mesh_.Sf()
    ),

    p(thermo_->p()),
    psi(thermo_->psi()),

    pressureReference_
    (
        thermo_->p(),
        pimple_.dict(),
        thermo_->incompressible()
    ),

    turbulence_
    (
        compressible::momentumTransportModel::New
        (
            rho,
            U,
            phi,
            thermo_()
        )
    ),

    thermophysicalTransport_
    (
        fluidThermophysicalTransportModel::New(turbulence_(), thermo_())
    ),

    dpdt
    (
        IOobject
        (
            "dpdt",
            runTime_.timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimPressure/dimTime, 0)
    ),

    K("K", 0.5*magSqr(U)),

    initialMass_(fvc::domainIntegrate(rho)),

    cumulativeContErr(0.0),

    MRF(mesh_),
    fvModels(fvModels::New(mesh_)),
    fvConstraints(fvConstraints::New(mesh_))
{
    {
        const Time& runTime = runTime_;
        #include "createRhoUfIfPresent.H"
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionSolvers::compressibleFluid::~compressibleFluid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::regionSolvers::compressibleFluid::solve()
{
    // --- Pressure-velocity PIMPLE corrector loop
    while (pimple_.loop())
    {
        // if (pimple_.firstPimpleIter())
        // {
        //     // Store momentum to set rhoUf for introduced faces.
        //     autoPtr<volVectorField> rhoU;
        //     if (rhoUf.valid())
        //     {
        //         rhoU = new volVectorField("rhoU", rho*U);
        //     }

        //     // Do any mesh changes
        //     mesh_.update();

        //     if (mesh_.changing())
        //     {
        //         MRF.update();

        //         // if (correctPhi)
        //         {
        //             // Calculate absolute flux
        //             // from the mapped surface velocity
        //             phi = mesh_.Sf() & rhoUf_();

        //             correctUphiBCs(rho, U, phi, true);

        //             CorrectPhi
        //             (
        //                 phi,
        //                 p,
        //                 rho,
        //                 psi,
        //                 dimensionedScalar("rAUf", dimTime, 1),
        //                 divrhoU_(),
        //                 pimple_
        //             );

        //             // Make the fluxes relative to the mesh-motion
        //             fvc::makeRelative(phi, rho, U);

        //         }

        //         // if (checkMeshCourantNo)
        //         // {
        //         //     #include "meshCourantNo.H"
        //         // }
        //     }
        // }

        if
        (
            !mesh_.steady()
         && !pimple_.simpleRho()
         && pimple_.firstPimpleIter()
        )
        {
            #include "rhoEqn.H"
        }

        // fvModels.correct();

        solveU();
        solveE();

        // --- Pressure corrector loop
        while (pimple_.correct())
        {
            solvep();
        }

        if (pimple_.turbCorr())
        {
            turbulence_->correct();
            thermophysicalTransport_->correct();
        }
    }

    if (!mesh_.steady())
    {
        rho = thermo_->rho();
    }
    return false;
}


void Foam::regionSolvers::compressibleFluid::solveU()
{
    MRF.correctBoundaryVelocity(U);

    tUEqn =
    (
        fvm::ddt(rho, U) + fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turbulence_->divDevTau(U)
     ==
        fvModels.source(rho, U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvConstraints.constrain(UEqn);

    if (pimple_.momentumPredictor())
    {
        ::Foam::solve(UEqn == -fvc::grad(p));

        fvConstraints.constrain(U);
        K = 0.5*magSqr(U);
    }
}


void Foam::regionSolvers::compressibleFluid::solveE()
{
    volScalarField& he = thermo_->he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div(fvc::absolute(phi, rho, U), p/rho)
          : -dpdt
        )
      + thermophysicalTransport_->divq(he)
     ==
        fvModels.source(rho, he)
    );

    EEqn.relax();

    fvConstraints.constrain(EEqn);

    EEqn.solve();

    fvConstraints.constrain(he);

    thermo_->correct();
}



void Foam::regionSolvers::compressibleFluid::solvep()
{
    if ((!mesh_.steady() && !pimple_.simpleRho()) || pimple_.consistent())
    {
        rho = thermo_->rho();
    }

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution
    const volScalarField psip0(psi*p);

    const volScalarField rAU("rAU", 1.0/tUEqn().A());
    const surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));

    tmp<volScalarField> rAtU
    (
        pimple_.consistent()
      ? volScalarField::New("rAtU", 1.0/(1.0/rAU - tUEqn().H1()))
      : tmp<volScalarField>(nullptr)
    );
    tmp<surfaceScalarField> rhorAtUf
    (
        pimple_.consistent()
      ? surfaceScalarField::New("rhoRAtUf", fvc::interpolate(rho*rAtU()))
      : tmp<surfaceScalarField>(nullptr)
    );

    const volScalarField& rAAtU = pimple_.consistent() ? rAtU() : rAU;
    const surfaceScalarField& rhorAAtUf =
        pimple_.consistent() ? rhorAtUf() : rhorAUf;

    volVectorField HbyA(constrainHbyA(rAU*tUEqn().H(), U, p));

    if (pimple_.nCorrPiso() <= 1)
    {
        tUEqn.clear();
    }

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::interpolate(rho)*fvc::flux(HbyA)
      + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi, rhoUf))
    );

    fvc::makeRelative(phiHbyA, rho, U);
    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    bool adjustMass = false;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, rho, U, phiHbyA, rhorAAtUf, MRF);

    if (pimple_.transonic())
    {
        surfaceScalarField phid
        (
            "phid",
            (fvc::interpolate(psi)/fvc::interpolate(rho))*phiHbyA
        );

        phiHbyA -=
            fvc::interpolate(psi*p)
           *phiHbyA
           /fvc::interpolate(rho);

        if (pimple_.consistent())
        {
            phiHbyA += (rhorAAtUf - rhorAUf)*fvc::snGrad(p)*mesh_.magSf();
            HbyA += (rAAtU - rAU)*fvc::grad(p);
        }

        fvScalarMatrix pDDtEqn
        (
            fvc::ddt(rho) + psi*correction(fvm::ddt(p))
          + fvc::div(phiHbyA) + fvm::div(phid, p)
         ==
            fvModels.source(psi, p, rho.name())
        );

        while (pimple_.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAAtUf, p));

            // Relax the pressure equation to ensure diagonal-dominance
            pEqn.relax();

            pEqn.setReference
            (
                pressureReference_.refCell(),
                pressureReference_.refValue()
            );

            pEqn.solve();

            if (pimple_.finalNonOrthogonalIter())
            {
                phi = phiHbyA + pEqn.flux();
            }
        }
    }
    else
    {
        if (mesh_.steady())
        {
            adjustMass = adjustPhi(phiHbyA, U, p);
        }

        if (pimple_.consistent())
        {
            phiHbyA += (rhorAAtUf - rhorAUf)*fvc::snGrad(p)*mesh_.magSf();
            HbyA += (rAAtU - rAU)*fvc::grad(p);
        }

        fvScalarMatrix pDDtEqn
        (
            fvc::ddt(rho) + psi*correction(fvm::ddt(p))
          + fvc::div(phiHbyA)
         ==
            fvModels.source(psi, p, rho.name())
        );

        while (pimple_.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAAtUf, p));

            pEqn.setReference
            (
                pressureReference_.refCell(),
                pressureReference_.refValue()
            );

            pEqn.solve();

            if (pimple_.finalNonOrthogonalIter())
            {
                phi = phiHbyA + pEqn.flux();
            }
        }
    }

    if (mesh_.steady())
    {
        const fvMesh& mesh = mesh_;
        const Time& runTime = runTime_;
        #include "continuityErrs.H"
    }
    else
    {
        const bool constrained = fvConstraints.constrain(p);

        // Thermodynamic density update
        thermo_->correctRho(psi*p - psip0);

        if (constrained)
        {
            rho = thermo_->rho();
        }

        #include "rhoEqn.H"

        const fluidThermo& thermo = thermo_;
        #include "compressibleContinuityErrs.H"
    }

    // Explicitly relax pressure for momentum corrector
    p.relax();

    U = HbyA - rAAtU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvConstraints.constrain(U);
    K = 0.5*magSqr(U);

    if (mesh_.steady())
    {
        fvConstraints.constrain(p);
    }

    // For steady compressible closed-volume cases adjust the pressure level
    // to obey overall mass continuity
    if (adjustMass && !thermo_->incompressible())
    {
        p += (initialMass_ - fvc::domainIntegrate(thermo_->rho()))
            /fvc::domainIntegrate(psi);
        p.correctBoundaryConditions();
    }

    if (mesh_.steady() || pimple_.simpleRho() || adjustMass)
    {
        rho = thermo_->rho();
    }

    // Correct rhoUf if the mesh is moving
    fvc::correctRhoUf(rhoUf, rho, U, phi);

    if ((mesh_.steady() || pimple_.simpleRho()) && !pimple_.transonic())
    {
        rho.relax();
    }

    if (thermo_->dpdt())
    {
        dpdt = fvc::ddt(p);

        if (mesh_.moving())
        {
            dpdt -= fvc::div(fvc::meshPhi(rho, U), p);
        }
    }
}

Foam::scalar Foam::regionSolvers::compressibleFluid::CoNum() const
{
    const fvMesh& mesh = mesh_;
    const Time& runTime = runTime_;
    #include "compressibleCourantNo.H"
    return CoNum;
}


Foam::scalar Foam::regionSolvers::compressibleFluid::maxCo() const
{
    return 1.0;
}

// ************************************************************************* //
