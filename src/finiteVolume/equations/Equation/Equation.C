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

#include "Equation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Equation<Type>::Equation()
:
    lowerLimit_(-great*pTraits<Type>::one),
    upperLimit_(great*pTraits<Type>::one)
{}


template<class Type>
Foam::Equation<Type>::Equation(const Type lowerLimit, const Type upperLimit)
:
    lowerLimit_(lowerLimit),
    upperLimit_(upperLimit)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Equation<Type>::~Equation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
bool Foam::Equation<Type>::checkBounds(const scalar x) const
{
    if (x < lowerLimit_ || x > upperLimit_)
    {
        #ifdef FULLDEBUG
        FatalErrorInFunction
            << "Request function evaluation is out of bounds." << nl
            << "lowerLimit: " << lowerLimit_ << endl
            << "upperLimit: " << upperLimit_ << endl
            << "x: " << x << endl
            << abort(FatalError);
        #endif
        return false;
    }
    return false;
}

// ************************************************************************* //
