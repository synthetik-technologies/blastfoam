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

#include "UnivariateEquation.H"

// * * * * * * * * * * * * * Static member functions * * * * * * * * * * * * //

// template<class Type>
// bool Foam::UnivariateEquation<Type>::checkJacobian() const
// {
//     if
//     (
//         (void*)(this->*(&UnivariateEquation<Type>::jacobian))
//      == (void*)(&UnivariateEquation<Type>::jacobian)
//     )
//     {
//         return false;
//     }
//     return true;
// }
//
//
// template<class Type>
// bool Foam::UnivariateEquation<Type>::checkHessian() const
// {
//     return false;
//     if
//     (
//         (void*)(this->*(&UnivariateEquation<Type>::hessian))
//      == (void*)(&UnivariateEquation<Type>::hessian)
//     )
//     {
//         return false;
//     }
//     return true;
// }


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationType>
Foam::UnivariateEquation<EquationType>::UnivariateEquation
(
    const label nVar,
    const inType& lowerLimit,
    const inType& upperLimit
)
:
    EquationType
    (
        nVar,
        1,
        lowerLimit,    
        upperLimit
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class EquationType>
Foam::UnivariateEquation<EquationType>::~UnivariateEquation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// template<class Type>
// Foam::label Foam::UnivariateEquation<Type>::nDerivatives() const
// {
//     label nDeriv = 0;
//
//     // Check if Jacobian has been implemented
//     if (checkJacobian())
//     {
//         nDeriv++;
//     }
//     else
//     {
//         return nDeriv;
//     }
//
//     Check if Hessian has been implemented
//     if (checkHessian())
//     {
//         nDeriv++;
//     }
//     return nDeriv;
// }


// ************************************************************************* //
