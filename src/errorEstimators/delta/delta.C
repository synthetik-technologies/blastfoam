/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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

#include "delta.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace errorEstimators
{
    defineTypeNameAndDebug(delta, 0);
    addToRunTimeSelectionTable(errorEstimator, delta, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::errorEstimators::delta::delta
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    errorEstimator(mesh, dict),
    name_(dict.lookup("deltaField"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::errorEstimators::delta::~delta()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::errorEstimators::delta::update()
{
    const volScalarField& x = mesh_.lookupObject<volScalarField>(name_);
    volScalarField& error(*this);

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const label nInternalFaces = mesh_.nInternalFaces();
    error = 0.0;

    for (label facei = 0; facei < nInternalFaces; facei++)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        scalar eT = mag(x[own] - x[nei])/Foam::min(x[own], x[nei]);
        error[own] = Foam::max(error[own], eT);
        error[nei] = Foam::max(error[nei], eT);
    }

    // Boundary faces
    forAll(error.boundaryField(), patchi)
    {
        if (error.boundaryField()[patchi].coupled())
        {
            const fvPatch& p = x.boundaryField()[patchi].patch();

            const labelUList& faceCells = p.faceCells();
            scalarField fp(x.boundaryField()[patchi].patchInternalField());

            scalarField fn
            (
                x.boundaryField()[patchi].patchNeighbourField()
            );

            forAll(faceCells, facei)
            {
                scalar eT =
                    mag(fp[facei] - fn[facei])/Foam::min(fp[facei], fn[facei]);
                error[faceCells[facei]]=
                    Foam::max(error[faceCells[facei]], eT);
            }
        }
    }
    error.correctBoundaryConditions();
}

// ************************************************************************* //
