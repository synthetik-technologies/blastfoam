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
    )
{
    this->read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::errorEstimators::volumeFraction::~volumeFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::errorEstimators::volumeFraction::read(const dictionary& dict)
{
    scalar threshold = dict.lookup<scalar>("threshold");
    lowerRefine_ = threshold;
    lowerUnrefine_ = threshold;
    upperRefine_ = 1.0 - threshold;
    upperUnrefine_ = 1.0 - threshold;

    if (dict.found("maxRefinement"))
    {
        maxLevel_ = dict.lookup<label>("maxRefinement");
        minDx_ = -1;
    }
    else if (dict.found("minDx"))
    {
        minDx_ = dict.lookup<scalar>("minDx");
        maxLevel_ = -1;
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Either maxRefinement or minDx must be specified" << endl
            << abort(FatalIOError);
    }
}

// ************************************************************************* //
