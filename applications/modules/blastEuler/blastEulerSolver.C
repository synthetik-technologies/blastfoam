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
    integrator_(mesh),
    g_
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedVector("g", dimAcceleration, vector::zero)
    ),
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
            surfaceScalarField::New
            (
                IOobject::groupName("amaxSf", phases[phasei].name()),
                mesh_,
                dimensionedScalar(dimVelocity*dimArea, Zero)
            )
        );

        tmp<volScalarField> tc(phases[phasei].speedOfSound());
        const volScalarField& c = tc();
        const volVectorField& U = phases[phasei].U();

        const scalarField& magSf = mesh_.magSf();
        const labelList& owner = mesh_.faceOwner();
        const labelList& neighbour = mesh_.faceNeighbour();
        forAll(neighbour, facei)
        {
            amaxSf[facei] =
                sqrt
                (
                    max
                    (
                        magSqr(U[owner[facei]]) + sqr(c[owner[facei]]),
                        magSqr(U[neighbour[facei]]) + sqr(c[neighbour[facei]])
                    )
                )*magSf[facei];
        }

        // Remove wave speed from wedge boundaries
        surfaceScalarField::Boundary& bamaxSf = amaxSf.boundaryFieldRef();
        forAll(amaxSf.boundaryField(), patchi)
        {
            const fvPatch& patch = mesh_.boundary()[patchi];
            const scalarField& pmagSf = patch.magSf();
            const labelList& faceCells = patch.faceCells();
            const fvPatchVectorField& pU = U.boundaryField()[patchi];
            const fvPatchScalarField& pc = c.boundaryField()[patchi];
            fvsPatchScalarField& pamaxSf = bamaxSf[patchi];
            if (patch.coupled())
            {
                const vectorField nbrU(pU.patchNeighbourField());
                const scalarField nbrc(pc.patchNeighbourField());
                forAll(pU, fi)
                {
                    const label own = faceCells[fi];
                    pamaxSf[fi] =
                        sqrt
                        (
                            max
                            (
                                magSqr(U[own]) + sqr(c[own]),
                                magSqr(nbrU[fi]) + sqr(nbrc[fi])
                            )
                        )*pmagSf[fi];
                }
            }
            else if (!isA<wedgeFvPatch>(patch) && !isA<emptyFvPatch>(patch))
            {
                if
                (
                    (mesh_.dynamic() || mesh_.distributing())
                 && (mesh_.moving())
                )
                {
                    forAll(pU, fi)
                    {
                        const label own = faceCells[fi];
                        pamaxSf[fi] = (mag(U[own]) + c[own])*pmagSf[fi];
                    }
                }
                else
                {
                    forAll(pU, fi)
                    {
                        pamaxSf[fi] = (mag(pU[fi]) + pc[fi])*pmagSf[fi];
                    }
                }
            }
        }

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
            0.5*gMax(sumAmaxSf/V)*mesh_.time().deltaTValue();
        meanPhaseCoNums[phasei] =
            0.5*(gSum(amaxSf)/gSum(V))*mesh_.time().deltaTValue();
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
    return max(phaseCoNums);
}


void Foam::solvers::blastEuler::solve()
{
    Info<< "Calculating Fluxes" << endl;
    integrator_.integrate();
}


void Foam::solvers::blastEuler::postSolve()
{
    fluid_.printInfo();
    integrator_.clear();
}


// ************************************************************************* //
