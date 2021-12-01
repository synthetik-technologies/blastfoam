/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020-2021
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

#include "multicomponent.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace errorEstimators
{
    defineTypeNameAndDebug(multicomponent, 0);
    addToRunTimeSelectionTable(errorEstimator, multicomponent, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::errorEstimators::multicomponent::multicomponent
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& name
)
:
    errorEstimator(mesh, dict, name),
    names_(),
    errors_()
{
    PtrList<entry> errorEntries(dict.lookup("errorEstimators"));
    errors_.resize(errorEntries.size());
    forAll(errorEntries, i)
    {
        names_.append(errorEntries[i].keyword());
        errors_.set
        (
            i,
            errorEstimator::New
            (
                mesh,
                errorEntries[i].dict(),
                errorEntries[i].keyword()
            ).ptr()
        );
    }
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::errorEstimators::multicomponent::~multicomponent()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::errorEstimators::multicomponent::read(const dictionary& dict)
{
    PtrList<entry> errorEntries(dict.lookup("errorEstimators"));
    forAll(errors_, i)
    {
        errors_[i].read(errorEntries[i].dict());
    }
    maxLevel_ = 0;
    forAll(names_, i)
    {
        maxLevel_ = max(maxLevel_, errors_[i].maxLevel());
    }
}


void Foam::errorEstimators::multicomponent::update(const bool scale)
{
    if (updateCurTimeIndex(!scale))
    {
        return;
    }

    volScalarField error
    (
        IOobject
        (
            "error",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        -1.0
    );

    forAll(errors_, i)
    {
        errors_[i].update(scale);
        error = max(error,  errors_[i].error());
    }
    error_ = error;
}


Foam::labelList Foam::errorEstimators::multicomponent::maxRefinement() const
{
    maxRefinement_.resize(mesh_.nCells());
    maxRefinement_ = 0;
    forAll(errors_, i)
    {
        const volScalarField& errori(errors_[i].error());
        labelList maxCellLevel(errors_[i].maxRefinement());
        forAll(maxRefinement_, celli)
        {
            if
            (
                errori[celli] > 0.5
             && maxRefinement_[celli] < maxCellLevel[celli]
            )
            {
                maxRefinement_[celli] = maxCellLevel[celli];
                const labelList& cellCells(mesh_.cellCells()[celli]);
                forAll(cellCells, j)
                {
                    label cellj = cellCells[j];
                    maxRefinement_[cellj] = maxCellLevel[celli];
                }

            }
        };
    }
    return maxRefinement_;
}

// ************************************************************************* //
