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

#include "ScalarEquation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class InType>
Foam::ScalarEquation<InType>::ScalarEquation
(
    const label nVar,
    const InType& lowerLimit,
    const InType& upperLimit
)
:
    Equation<InType, scalar>(nVar, 1, lowerLimit, upperLimit)
{}


template<class InType>
Foam::ScalarEquation<InType>::ScalarEquation
(
    const label nVar,
    const label,
    const InType& lowerLimit,
    const InType& upperLimit
)
:
    Equation<InType, scalar>(nVar, 1, lowerLimit, upperLimit)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class InType>
Foam::ScalarEquation<InType>::~ScalarEquation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class InType>
void Foam::ScalarEquation<InType>::calculateGradient
(
    const InType& x0,
    const label li,
    const scalar& fx0,
    scalarField& grad
) const
{
    scalar fx1;
    scalarField x1(x0);
    for (label cmpti = 0; cmpti < this->nVar_; cmpti++)
    {
        x1 = x0;
        x1[cmpti] += this->dx(cmpti);
        this->f(x1, li, fx1);
        grad[cmpti] = (fx1 - fx0)/this->dx(cmpti);
    }
}


template<class InType>
void Foam::ScalarEquation<InType>::gradient
(
    const InType& x0,
    const label li,
    scalar& fx0,
    scalarField& grad
) const
{
    Equation<InType, scalar>::f(x0, li, fx0);
    calculateGradient(x0, li, fx0, grad);
}


// ************************************************************************* //
