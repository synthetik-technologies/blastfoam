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

#include "minimizationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::minimizationScheme>
Foam::minimizationScheme::New
(
    const word& minimizationSchemeType,
    const scalarUnivariateEquation& eqn,
    const dictionary& dict
)
{
    Info<< "Selecting minimization scheme: " << minimizationSchemeType
        << endl;
    const dictionary& coeffDict = dict.optionalSubDict
    (
        minimizationSchemeType
    );
    if
    (
        isA<scalarEquation>(eqn)
      && !dictionaryMultivariateConstructorTablePtr_->found(minimizationSchemeType)
    )
    {

        dictionaryUnivariateConstructorTable::iterator cstrIter =
            dictionaryUnivariateConstructorTablePtr_->find(minimizationSchemeType);

        if (cstrIter == dictionaryUnivariateConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown univariate minimization scheme type "
                << minimizationSchemeType << nl << nl
                << "Valid univariate minimization schemes are : " << endl
                << dictionaryUnivariateConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }
        return autoPtr<minimizationScheme>
        (
            cstrIter()(eqn, coeffDict)
        );
    }

    dictionaryMultivariateConstructorTable::iterator cstrIter =
        dictionaryMultivariateConstructorTablePtr_->find(minimizationSchemeType);

    if (cstrIter == dictionaryMultivariateConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown multivariate minimization scheme type "
            << minimizationSchemeType << nl << nl
            << "Valid multivariate minimization schemes:" << endl
            << dictionaryMultivariateConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
    return autoPtr<minimizationScheme>
    (
        cstrIter()(eqn, coeffDict)
    );
}

Foam::autoPtr<Foam::minimizationScheme>
Foam::minimizationScheme::New
(
    const scalarUnivariateEquation& eqn,
    const dictionary& dict
)
{
    return New(dict.lookup<word>("solver"), eqn, dict);
}
// ************************************************************************* //
