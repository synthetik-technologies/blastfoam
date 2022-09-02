/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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
    RBFFunction

\*---------------------------------------------------------------------------*/

#include "RBFFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<RBFFunction> RBFFunction::New
(
    const dictionary& dict
)
{
    const word type(dict.lookup("RBFFunction"));
    Info<< "Radial Basis Function interpolation: "
        << "Selecting RBF function: " << type << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown RBFFunction type "
            << type << endl << endl
            << "Valid  RBFFunctions are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << abort(FatalIOError);
    }

    return autoPtr<RBFFunction>
    (
        cstrIter()(dict.optionalSubDict(type + "Coeffs"))
    );
}

autoPtr<RBFFunction> RBFFunction::New
(
    const word& type,
    const dictionary& dict
)
{
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown RBFFunction type "
            << type << endl << endl
            << "Valid  RBFFunctions are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << abort(FatalIOError);
    }

    return autoPtr<RBFFunction>
    (
        cstrIter()(dict.optionalSubDict(type + "Coeffs"))
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
