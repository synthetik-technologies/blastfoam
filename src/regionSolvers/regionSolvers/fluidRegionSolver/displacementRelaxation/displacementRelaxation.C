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

#include "displacementRelaxation.H"
#include "valuePointPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(displacementRelaxation, 0);
    defineRunTimeSelectionTable(displacementRelaxation, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementRelaxation::displacementRelaxation
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    coupledPatches_(globalPolyBoundaryMesh::New(mesh_).coupledPatches()),
    initialError_(-1),
    error_(great)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacementRelaxation::~displacementRelaxation()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::displacementRelaxation::updateError
(
    const pointVectorField& p
)
{
    error_ = Zero;
    labelList coupledPatches(globalPolyBoundaryMesh::New(mesh_).coupledPatches());
    forAll(coupledPatches, pi)
    {
        const label patchi = coupledPatches[pi];
        const pointPatchVectorField& pp = p.boundaryField()[patchi];
        const pointPatchVectorField& ppPrev =
            initialError_ < 0
          ? p.oldTime().boundaryField()[patchi]
          : p.prevIter().boundaryField()[patchi];
        if (isA<valuePointPatchVectorField>(pp))
        {
            error_ +=
                sum
                (
                    magSqr
                    (
                        dynamicCast<const valuePointPatchVectorField>(pp)
                      - dynamicCast<const valuePointPatchVectorField>(ppPrev)
                    )
                );
        }
    }
    reduce(error_, sumOp<scalar>());
    error_ = sqrt(error_);

    if (initialError_ < 0)
    {
        initialError_ = error_;
    }
}


void Foam::displacementRelaxation::clear()
{
    initialError_ = -1;
    error_ = great;
}


// ************************************************************************* //
