/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    crackPathLimiter

\*---------------------------------------------------------------------------*/

#include "noCrackPathLimiter.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<crackPathLimiter> crackPathLimiter::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    word lawTypeName(crackPathLimiters::none::typeName);
    const dictionary* coeffDictPtr = &dict;
    if (dict.isDict(typeName))
    {
        coeffDictPtr = &dict.subDict(typeName);
        lawTypeName = coeffDictPtr->lookup<word>("type");
    }
    else if (dict.found(typeName))
    {
        lawTypeName = dict.lookup<word>(typeName);
        coeffDictPtr = &dict.optionalSubDict(lawTypeName + "Coeffs");
    }

    Info<< "Selecting crack path limiter: " << lawTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(lawTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown crackPathLimiter type "
            << lawTypeName << endl << endl
            << "Valid  crackPathLimiters are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return autoPtr<crackPathLimiter>(cstrIter()(mesh, *coeffDictPtr));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
