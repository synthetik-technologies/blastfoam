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

#include "fluidRegionSolver.H"
#include "surfaceFields.H"
#include "fvc.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionSolvers
{
    defineTypeNameAndDebug(fluid, 0);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionSolvers::fluid::fluid
(
    fvMesh& mesh,
    const regionSolverList& regions
)
:
    regionSolver(mesh, regions),
    velocityFields_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionSolvers::fluid::~fluid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionSolvers::fluid::initialiseMesh(const IterType iter)
{
    moveMesh(FINAL_ITER);

    accelerationSchemes_.clear();

    if (mesh_.moving())
    {
        const_cast<surfaceScalarField&>(mesh_.phi()) == Zero;
        forAll(velocityFields_, i)
        {
            if (mesh_.foundObject<volVectorField>(velocityFields_[i]))
            {
                mesh_.lookupObjectRef<volVectorField>
                (
                    velocityFields_[i]
                ).correctBoundaryConditions();
            }
        }
    }

    if (mesh_.moving() && iter == FINAL_ITER)
    {
        mesh_.lookupObject<pointIOField>("points").write();
    }
}


void Foam::regionSolvers::fluid::initialise()
{
    moveMesh(FINAL_ITER);

    // Set old points
    if (mesh_.moving())
    {
        const_cast<surfaceScalarField&>(mesh_.phi()) == Zero;
        forAll(velocityFields_, i)
        {
            if (mesh_.foundObject<volVectorField>(velocityFields_[i]))
            {
                mesh_.lookupObjectRef<volVectorField>
                (
                    velocityFields_[i]
                ).correctBoundaryConditions();
            }
        }
        mesh_.resetMotion();
    }
}


bool Foam::regionSolvers::fluid::changeMesh()
{
    bool changed = regionSolver::changeMesh();
    return changed;
}


bool Foam::regionSolvers::fluid::moveMesh(const IterType iter)
{
    regionSolver::moveMesh(iter);

    if (mesh_.moving())
    {
        storePrevIter();

        if (mesh_.moving() && (debug || regionSolver::debug))
        {
            Info<<"Mesh boundary velocity (max/mean): " << endl;
            forAll(mesh_.boundary(), patchi)
            {
                const fvPatch& p = mesh_.boundary()[patchi];
                if (!p.coupled() && returnReduce(p.size(), sumOp<label>()))
                {
                    const polyPatch& pp = p.patch();
                    const pointField& oldPoints = mesh_.oldPoints();

                    vectorField oldFc(pp.size());
                    forAll(oldFc, i)
                    {
                        oldFc[i] = pp[i].centre(oldPoints);
                    }

                    const scalar deltaT = mesh_.time().deltaTValue();

                    const vectorField Up((pp.faceCentres() - oldFc)/deltaT);

                    const volVectorField& U =
                        mesh_.lookupObject<volVectorField>("U");

                    scalarField phip
                    (
                        p.patchField<surfaceScalarField, scalar>(fvc::meshPhi(U))
                    );

                    const vectorField n(p.nf());
                    const scalarField& magSf = p.magSf();
                    scalarField Un(phip/(magSf + vSmall));
                    const vectorField pU(Up + n*(Un - (n & Up)));

                    Info<< "    " << mesh_.boundary()[patchi].name()<<": "
                        << gMaxMagSqr(pU) << "/" << gAverage(pU) << endl;
                }
            }
        }
    }

    accelerationSchemes_.print(Info);

    if (mesh_.moving())
    {
        forAll(velocityFields_, i)
        {
            if (mesh_.foundObject<volVectorField>(velocityFields_[i]))
            {
                mesh_.lookupObjectRef<volVectorField>
                (
                    velocityFields_[i]
                ).correctBoundaryConditions();
            }
        }
    }
    return mesh_.moving();
}


// ************************************************************************* //
