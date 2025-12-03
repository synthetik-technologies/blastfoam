/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2025
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "blastEulerSolver.H"
#include "wedgeFvPatch.H"
#include "emptyFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(blastEuler, 0);
    addToRunTimeSelectionTable(solver, blastEuler, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::blastEuler::blastEuler(fvMesh& mesh)
:
    explicitSolver(mesh),
    integrator_(mesh, false),
    fluid_(mesh)
{
    integrator_.addSystem(fluid_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::blastEuler::~blastEuler()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::solvers::blastEuler::CoNum() const
{
    const PtrList<phaseModel>& phases = fluid_.phases();
    const scalarField& V = mesh_.V();

    scalarField sumPhi(V.size(), 0.0);
    scalarList phaseCoNums(phases.size(), 0.0);
    scalarList meanPhaseCoNums(phases.size(), 0.0);
    forAll(phases, phasei)
    {
        surfaceScalarField amaxSf
        (
            fvc::interpolate(phases[phasei].speedOfSound())*mesh_.magSf()
        );

        // Remove wave speed from wedge boundaries
        forAll(amaxSf.boundaryField(), patchi)
        {
            if (isA<wedgeFvPatch>(mesh_.boundary()[patchi]))
            {
                amaxSf.boundaryFieldRef() = Zero;
            }
        }
        amaxSf += mag(fvc::flux(phases[phasei].U()));

        scalarField sumAmaxSf
        (
            fvc::surfaceSum(amaxSf)().primitiveField()
        );

        sumPhi +=
            fvc::surfaceSum
            (
                amaxSf*fvc::interpolate(phases[phasei])
            )().primitiveField();

        phaseCoNums[phasei] =
            0.5*gMax(sumAmaxSf/V)*runTime.deltaTValue();
        meanPhaseCoNums[phasei] =
            0.5*(gSum(amaxSf)/gSum(V))*runTime.deltaTValue();
    }

    scalar maxCoNum = 0.5*gMax(sumPhi/V)*mesh_.time().deltaTValue();
    scalar meanCoNum
    (
        0.5*(gSum(sumPhi)/gSum(V))*mesh_.time().deltaTValue()
    );

    Info<< "Courant number: mean = " << meanCoNum
        << ", max = " << maxCoNum << endl;

    Info<< "Phase Courant numbers based on eigenvalues:"
        << incrIndent << endl;

    forAll(phases, phasei)
    {
        Info<< indent << phases[phasei].name() << ": "
            << "mean = " << meanPhaseCoNums[phasei]
            << ", max = " << phaseCoNums[phasei] << nl;
    }
    Info<< endl << decrIndent;

    // bool hasMassTransfer = false;
    // forAll(phases, phasei)
    // {
    //     const phaseModel& phase = phases[phasei];
    //
    //     if (fluid.hasMassTransfer(phase))
    //     {
    //         hasMassTransfer = true;
    //         break;
    //     }
    // }

    // mDotCoNum = 0.0;
    // if (hasMassTransfer)
    // {
    //     Info<< "Maximum phase Courant numbers based on mass transfer:" << endl
    //         << incrIndent;
    //
    //     const dimensionedScalar zeroMDot(dimDensity/dimTime, 0.0);
    //     forAll(phases, phasei)
    //     {
    //         const phaseModel& phase = phases[phasei];
    //
    //         if (fluid.hasMassTransfer(phase))
    //         {
    //             volScalarField totalMDot
    //             (
    //                 volScalarField::New
    //                 (
    //                     IOobject::groupName("mDot", phase.name()),
    //                     mesh,
    //                     dimensionedScalar(dimDensity/dimTime, 0.0)
    //                 )
    //             );
    //             forAll(phases, phasej)
    //             {
    //                 const phaseModel& otherPhase = phases[phasej];
    //                 if (&otherPhase != &phase)
    //                 {
    //                     totalMDot += fluid.mDot(phase, otherPhase);
    //                 }
    //             }
    //             scalar mDotCo =
    //                 mag
    //                 (
    //                     gMaxMagSqr
    //                     (
    //                         (
    //                             totalMDot*runTime.deltaTValue()
    //                            /max(phase.alphaRho(), phase.residualAlphaRho())
    //                         )()
    //                     )
    //                 );
    //             Info<< indent << phase.name() << ": " << mDotCo << nl;
    //             mDotCoNum = max(mDotCoNum, mDotCo);
    //         }
    //     }
    //     Info<< endl << decrIndent;
    // }
    return max(phaseCoNums);
}


Foam::scalar Foam::solvers::blastEuler::DiNum() const
{
    return 0.0;
}


void Foam::solvers::blastEuler::solveExplicit()
{
    Info<< "Calculating Fluxes" << endl;
    integrator_.integrate
    (
        true,   // doExplicit
        true,   // doStore
        false,  // doImplicit
        false,  // doPost
        false   // doClear
    );
}


void Foam::solvers::blastEuler::solveImplicit()
{
    integrator_.solveImplicit();
}


void Foam::solvers::blastEuler::postSolve()
{
    integrator_.postUpdate();
    integrator_.clear();

    fluid_.printInfo();
}


// ************************************************************************* //
