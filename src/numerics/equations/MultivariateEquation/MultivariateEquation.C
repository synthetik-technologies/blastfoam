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
void Foam::MultivariateEquation<Type>::checkLimits() const
{
    if (lowerLimits_.size() != nEqns() || upperLimits_.size() != nEqns())
    {
        FatalErrorInFunction
            << "Limits have not been set, but are required for the " << nl
            << "requested root solver." << endl
            << abort(FatalError);
    }
}


template<class Type>
void Foam::MultivariateEquation<Type>::calculateJacobian
(
    const scalarList& x0,
    const label li,
    const List<Type>& f0,
    RectangularMatrix<Type>& J
) const
{
    J.setSize(nEqns(), x0.size());
    forAll(x0, cmptj)
    {
        scalarList x1(x0);
        x1[cmptj] += dx_[cmptj];

        scalarField f1(this->f(x1, li));

        for (label cmpti = 0; cmpti < nEqns(); cmpti++)
        {
            J(cmpti, cmptj) =
                (f1[cmpti] - f0[cmpti])/(x1[cmptj] - x0[cmptj]);
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::MultivariateEquation<Type>::MultivariateEquation(const label n)
:
    lowerLimits_(n, -great),
    upperLimits_(n, great),
    dx_(n, 1e-6)
{}


template<class Type>
Foam::MultivariateEquation<Type>::MultivariateEquation
(
    const scalarList& lowerLimits,
    const scalarList& upperLimits
)
:
    lowerLimits_(lowerLimits),
    upperLimits_(upperLimits),
    dx_(lowerLimits.size(), 1e-6)
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
bool Foam::MultivariateEquation<Type>::checkBounds(const scalarList& xs) const
{
    checkLimits();
    forAll(xs, i)
    {
        if (xs[i] < lowerLimits_[i] || xs[i] > upperLimits_[i])
        {
            #ifdef FULLDEBUG
            FatalErrorInFunction
                << "Request function evaluation is out of bounds." << nl
                << "lowerLimits: " << lowerLimits_ << endl
                << "upperLimits: " << upperLimits_ << endl
                << "x: " << xs << endl
                << abort(FatalError);
            #endif
            return false;
        }
    }
    return true;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::MultivariateEquation<Type>::f
(
    const scalarList& x,
    const label li
) const
{
    tmp<Field<Type>> tmpFx(new Field<Type>(this->nEqns()));
    this->f(x, li, tmpFx.ref());
    return tmpFx;
}


template<class Type>
void Foam::MultivariateEquation<Type>::jacobian
(
    const scalarList& x,
    const label li,
    List<Type>& fx,
    RectangularMatrix<Type>& J
) const
{
    fx.resize(this->nEqns());
    this->f(x, li, fx);
    calculateJacobian(x, li, fx, J);
}


// ************************************************************************* //
