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

#include "univariateMinimizationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::univariateMinimizationScheme> Foam::univariateMinimizationScheme::New
(
    const scalarUnivariateEquation& eqn,
    const dictionary& dict
)
{
    word univariateMinimizationSchemeTypeName(dict.lookup("solver"));
    label nDeriv = eqn.nDerivatives();
    DebugInfo
        << "Selecting root solver "
        << univariateMinimizationSchemeTypeName << endl;

    if (nDeriv <= 0)
    {
        dictionaryZeroConstructorTable::iterator cstrIter =
            dictionaryZeroConstructorTablePtr_->find
            (
                univariateMinimizationSchemeTypeName
            );

        if (cstrIter == dictionaryZeroConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown univariateMinimizationScheme type "
                << univariateMinimizationSchemeTypeName << nl << nl
                << "Valid univariateMinimizationSchemes for no derivatives are : " << endl
                << dictionaryZeroConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }
        return autoPtr<univariateMinimizationScheme>(cstrIter()(eqn, dict));
    }
    else if (nDeriv == 1)
    {
        dictionaryOneConstructorTable::iterator cstrIter =
            dictionaryOneConstructorTablePtr_->find(univariateMinimizationSchemeTypeName);

        if (cstrIter == dictionaryOneConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown univariateMinimizationScheme type "
                << univariateMinimizationSchemeTypeName << nl << nl
                << "Valid univariateMinimizationSchemes for one derivative are : " << endl
                << dictionaryOneConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }
        return autoPtr<univariateMinimizationScheme>(cstrIter()(eqn, dict));
    }
    dictionaryTwoConstructorTable::iterator cstrIter =
        dictionaryTwoConstructorTablePtr_->find(univariateMinimizationSchemeTypeName);

    if (cstrIter == dictionaryTwoConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown univariateMinimizationScheme type "
            << univariateMinimizationSchemeTypeName << nl << nl
            << "Valid univariateMinimizationSchemes for are : " << endl
            << dictionaryTwoConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
    return autoPtr<univariateMinimizationScheme>(cstrIter()(eqn, dict));
}


// ************************************************************************* //
