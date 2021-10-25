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

#include "MultivariateEquation.H"

// * * * * * * * * * * * * * Static member functions * * * * * * * * * * * * //

// template<class Type>
// bool Foam::MultivariateEquation<Type>::checkJacobian() const
// {
//     if
//     (
//         (void*)(this->*(&MultivariateEquation<Type>::jacobian))
//      == (void*)(&MultivariateEquation<Type>::jacobian)
//     )
//     {
//         return false;
//     }
//     return true;
// }
//
//
// template<class Type>
// bool Foam::MultivariateEquation<Type>::checkHessian() const
// {
//     return false;
//     if
//     (
//         (void*)(this->*(&MultivariateEquation<Type>::hessian))
//      == (void*)(&MultivariateEquation<Type>::hessian)
//     )
//     {
//         return false;
//     }
//     return true;
// }


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::MultivariateEquation<Type>::calculateJacobian
(
    const scalarField& x0,
    const label li,
    const Field<Type>& f0,
    RectangularMatrix<Type>& J
) const
{
    J.setSize(this->nEqns_, this->nVar_);
    for (label cmptj = 0; cmptj < this->nVar_; cmptj++)
    {
        scalarField x1(x0);
        x1[cmptj] += this->dx(cmptj);
        scalarField f1(this->f(x1, li));
        for (label cmpti = 0; cmpti < this->nEqns_; cmpti++)
        {
            J(cmpti, cmptj) =
                (f1[cmpti] - f0[cmpti])/(x1[cmptj] - x0[cmptj]);
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::MultivariateEquation<Type>::MultivariateEquation
(
    const label nEqns,
    const scalarField& lowerLimits,
    const scalarField& upperLimits
)
:
    Equation<scalarField, Field<Type>>
    (
        lowerLimits.size(),
        nEqns,
        lowerLimits,    
        upperLimits
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::MultivariateEquation<Type>::~MultivariateEquation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// template<class Type>
// Foam::label Foam::MultivariateEquation<Type>::nDerivatives() const
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

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::MultivariateEquation<Type>::f
(
    const scalarField& x,
    const label li
) const
{
    tmp<Field<Type>> tmpFx(new Field<Type>(this->nEqns_));
    this->f(x, li, tmpFx.ref());
    return tmpFx;
}


template<class Type>
void Foam::MultivariateEquation<Type>::jacobian
(
    const scalarField& x,
    const label li,
    Field<Type>& fx,
    RectangularMatrix<Type>& J
) const
{
    fx.resize(this->nEqns_);
    this->f(x, li, fx);
    calculateJacobian(x, li, fx, J);
}


// ************************************************************************* //
