/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
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

#include "lineSearch.H"
#include "approximateLineSearch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lineSearch, 0);
    defineRunTimeSelectionTable(lineSearch, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lineSearch::lineSearch
(
    const scalarUnivariateEquation& eqns,
    const dictionary& dict
)
:
    eqns_(eqns),
    lsEqn_
    (
        eqns,
        dict.lookupOrDefault<word>
        (
            "descentMethod",
            "PolakRibiere"
        )
    )
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::lineSearch::~lineSearch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::lineSearch>
Foam::lineSearch::New
(
    const scalarUnivariateEquation& eqn,
    const dictionary& dict
)
{
    const word lineSearchType(dict.lookupOrDefault<word>("lineSearch", "exact"));
    Info<< "Selecting lineSearch: " << lineSearchType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(lineSearchType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown lineSearch type "
            << lineSearchType << nl << nl
            << "Valid lineSearches are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
    return autoPtr<lineSearch>(cstrIter()(eqn, dict));
}


// ************************************************************************* //
