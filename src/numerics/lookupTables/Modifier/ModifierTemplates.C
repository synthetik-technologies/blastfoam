/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
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

#include "Modifier.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::Modifier<Type>> Foam::Modifier<Type>::New(const word& mod)
{
    if (mod == "none")
    {
        return autoPtr<Modifier<Type>>(new modifiers::None<Type>());
    }
    else
    {
        FatalErrorInFunction
            << mod << " is not a valid mod scheme" << nl
            << "Options are: " << nl
            << "    none" << nl
            << abort(FatalError);
    }
    return autoPtr<Modifier<Type>>();
}

// ************************************************************************* //
