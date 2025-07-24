/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2023
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

#include "volumeFractionErrorEstimator.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace errorEstimators
{
    defineTypeNameAndDebug(volumeFraction, 0);
    addToRunTimeSelectionTable(errorEstimator, volumeFraction, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::errorEstimators::volumeFraction::volumeFraction
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& name
)
:
    fieldValue
    (
        mesh,
        dict,
        name,
        IOobject::groupName
        (
            "alpha",
            dict.lookup<word>("phase")
        )
    ),
    refineFaces_(dict.lookupOrDefault("refineFaces", true))
{
    this->read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::errorEstimators::volumeFraction::~volumeFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::errorEstimators::volumeFraction::update(const bool scale)
{
    if (updateCurTimeIndex(!scale))
    {
        return;
    }

    const labelHashSet& eCells = this->errorCells();

    const volScalarField& alpha =
        mesh_.lookupObjectRef<volScalarField>(fieldName_);

    // Error is volume fraction
    error_ = alpha;

    // Average of lower and upper refinement values
    const scalar avgRefine = 0.5*(lowerRefine_ + upperRefine_);

    if (!refineFaces_)
    {
         if (scale)
        {
            normalize(error_, eCells);
        }
        return;
    }

    // Check own/nei of faces for an interface
    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        const label own = mesh_.faceOwner()[facei];
        const label nei = mesh_.faceNeighbour()[facei];
        const bool foundOwn = eCells.found(own);
        const bool foundNei = eCells.found(nei);

        if
        (
            (foundOwn || foundNei)
         && (
                (alpha[own] < lowerRefine_ && alpha[nei] > upperRefine_)
             || (alpha[own] > upperRefine_ && alpha[nei] < lowerRefine_)
            )
        )
        {
            if (foundOwn) error_[own] = avgRefine;
            if (foundNei) error_[nei] = avgRefine;
        }
    }

    volScalarField::Boundary& berror = error_.boundaryFieldRef();
    forAll(berror, patchi)
    {
        fvPatchScalarField& perror = berror[patchi];
        if (perror.coupled())
        {
            const fvPatchScalarField& palpha = alpha.boundaryField()[patchi];
            const labelList& faceCells = perror.patch().faceCells();
            const scalarField alphaNbr(palpha.patchNeighbourField());
            forAll(perror, fi)
            {
                if
                (
                    eCells.found(faceCells[fi])
                 && (
                        (
                            palpha[fi] < lowerRefine_
                         && alphaNbr[fi] > upperRefine_
                        )
                     || (
                            palpha[fi] > upperRefine_
                         && alphaNbr[fi] < lowerRefine_
                        )
                    )
                )
                {
                    error_[faceCells[fi]] = avgRefine;
                    perror[fi] = avgRefine;
                }
            }
        }
    }

    if (scale)
    {
        normalize(error_, eCells);
    }
}


void Foam::errorEstimators::volumeFraction::read(const dictionary& dict)
{
    scalar threshold = dict.lookup<scalar>("threshold");
    lowerRefine_ = threshold;
    lowerUnrefine_ = threshold;
    upperRefine_ = 1.0 - threshold;
    upperUnrefine_ = 1.0 - threshold;

    readCellZones(dict);
    readMaxRefinement(dict);
}

// ************************************************************************* //
