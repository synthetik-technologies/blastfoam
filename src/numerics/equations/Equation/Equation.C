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

// * * * * * * * * * * * * * Private member functions  * * * * * * * * * * * //

// template<class Type>
// bool Foam::Equation<Type>::checkFirstDerivative() const
// {
//     if
//     (
//         (void*)(this->*(&Equation<Type>::dfdx))
//      == (void*)(&Equation<Type>::dfdx)
//     )
//     {
//         return false;
//     }
//     return true;
// }
//
//
// template<class Type>
// bool Foam::Equation<Type>::checkSecondDerivative() const
// {
//     if
//     (
//         (void*)(this->*(&Equation<Type>::d2fdx2))
//      == (void*)(&Equation<Type>::d2fdx2)
//     )
//     {
//         return false;
//     }
//     return true;
// }
//
//
// template<class Type>
// bool Foam::Equation<Type>::checkThirdDerivative() const
// {
//     if
//     (
//         (void*)(this->*(&Equation<Type>::d3fdx3))
//      == (void*)(&Equation<Type>::d3fdx3)
//     )
//     {
//         return false;
//     }
//     return true;
// }
//
//
// template<class Type>
// bool Foam::Equation<Type>::checkFourthDerivative() const
// {
//     if
//     (
//         (void*)(this->*(&Equation<Type>::d4fdx4))
//      == (void*)(&Equation<Type>::d4fdx4)
//     )
//     {
//         return false;
//     }
//     return true;
// }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Equation<Type>::Equation()
:
    lowerLimit_(-great),
    upperLimit_(great)
{}


template<class Type>
Foam::Equation<Type>::Equation
(
    const scalar lowerLimit,
    const scalar upperLimit
)
:
    lowerLimit_(lowerLimit),
    upperLimit_(upperLimit)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Equation<Type>::~Equation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// template<class Type>
// Foam::label Foam::Equation<Type>::nDerivatives() const
// {
//     label nDeriv = 0;
//     // Check if first derivative has been implemented
//     if (checkFirstDerivative())
//     {
//         nDeriv++;
//     }
//     else
//     {
//         return nDeriv;
//     }
//
//     // Check if second derivative has been implemented
//     if (checkSecondDerivative())
//     {
//         nDeriv++;
//     }
//     else
//     {
//         return nDeriv;
//     }
//
//     // Check if third derivative has been implemented
//     if (checkThirdDerivative())
//     {
//         nDeriv++;
//     }
//     else
//     {
//         return nDeriv;
//     }
//
//     // Check if fourth derivative has been implemented
//     if (checkFourthDerivative())
//     {
//         nDeriv++;
//     }
//     return nDeriv;
// }


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
