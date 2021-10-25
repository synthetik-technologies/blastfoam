/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
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

#include "minimizationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::minimizationScheme>
Foam::minimizationScheme::New
(
    const scalarEquation& eqns,
    const dictionary& dict
)
{
    word minimizationSchemeTypeName(dict.lookup("solver"));
    Info<< "Selecting root solver " << minimizationSchemeTypeName << endl;

    if (eqns.nVar() == 1)
    {
        dictionaryUnivariateConstructorTable::iterator cstrIter =
            dictionaryUnivariateConstructorTablePtr_->find(minimizationSchemeTypeName);

        if (cstrIter == dictionaryUnivariateConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown minimizationScheme type "
                << minimizationSchemeTypeName << nl << nl
                << "Valid minimizationSchemes for no derivatives are : " << endl
                << dictionaryUnivariateConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }
        return autoPtr<minimizationScheme>(cstrIter()(eqns, dict));
    }

    dictionaryMultivariateConstructorTable::iterator cstrIter =
        dictionaryMultivariateConstructorTablePtr_->find(minimizationSchemeTypeName);

    if (cstrIter == dictionaryMultivariateConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown minimizationScheme type "
            << minimizationSchemeTypeName << nl << nl
            << "Valid minimizationSchemes for one derivative are : " << endl
            << dictionaryMultivariateConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
    return autoPtr<minimizationScheme>(cstrIter()(eqns, dict));
}


// ************************************************************************* //
