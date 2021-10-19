/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2020
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

#include "scaledDelta.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace errorEstimators
{
    defineTypeNameAndDebug(scaledDelta, 0);
    addToRunTimeSelectionTable(errorEstimator, scaledDelta, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::errorEstimators::scaledDelta::scaledDelta
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& name
)
:
    errorEstimator(mesh, dict, name),
    fieldName_(dict.lookup("scaledDeltaField")),
    minVal_(dict.lookupOrDefault<scalar>("minValue", small))
{
    this->read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::errorEstimators::scaledDelta::~scaledDelta()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::errorEstimators::scaledDelta::update(const bool scale)
{
    volScalarField x
    (
        IOobject
        (
            "mag(" + fieldName_ + ")",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        0.0
    );
    this->getFieldValue(fieldName_, x);

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const label nInternalFaces = mesh_.nInternalFaces();
    error_ = 0.0;

    for (label facei = 0; facei < nInternalFaces; facei++)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        scalar eT =
            mag(x[own] - x[nei])/max(min(x[own], x[nei]), minVal_);
        error_[own] = max(error_[own], eT);
        error_[nei] = max(error_[nei], eT);
    }

    // Boundary faces
    forAll(error_.boundaryField(), patchi)
    {
        if (error_.boundaryField()[patchi].coupled())
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
                    mag(fp[facei] - fn[facei])
                   /max(min(fp[facei], fn[facei]), minVal_);
                error_[faceCells[facei]] =
                    max(error_[faceCells[facei]], eT);
            }
        }
    }

    if (scale)
    {
        normalize(error_);
    }
}

// ************************************************************************* //
