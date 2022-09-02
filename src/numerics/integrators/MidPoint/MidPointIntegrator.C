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

#include "MidPointIntegrator.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::MidPointIntegrator<Type>::MidPointIntegrator
(
    const equationType& eqn,
    const dictionary& dict
)
:
    Integrator<Type>(eqn, dict)
{}


template<class Type>
Foam::MidPointIntegrator<Type>::MidPointIntegrator
(
    const equationType& eqn,
    const integrator& inter
)
:
    Integrator<Type>(eqn, inter)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::MidPointIntegrator<Type>::integrate_
(
    const Type& Q,
    const scalar x0,
    const scalar x1,
    const scalar tol,
    const label li
) const
{
    const scalar dx(x1 - x0);
    if (mag(dx) <= this->minDx_)
    {
        return Q;
    }

    this->intervals_++;
    const scalar xm = 0.5*(x1 + x0);
    const scalar x0m = 0.5*(x0 + xm);
    const Type fx0((xm - x0)*this->eqnPtr_->fx(x0m, li));

    const scalar xm1 = 0.5*(xm + x1);
    const Type fx1((x1 - xm)*this->eqnPtr_->fx(xm1, li));
    this->evals_ += 2;

    const Type fx(fx0 + fx1);
    if (this->converged(fx, Q, dx, tol))
    {
        return fx;
    }
    else
    {
        return
            integrate_(fx0, x0, xm, tol/2.0, li)
          + integrate_(fx1, xm, x1, tol/2.0, li);
    }
}


template<class Type>
Type Foam::MidPointIntegrator<Type>::integrate
(
    const scalar X0,
    const scalar X1,
    const label li
) const
{
    scalar dx(X1 - X0);
    if (mag(dx) < small)
    {
        return dx*this->eqnPtr_->fx(X0, li);
    }

    this->reset(dx);

    if (this->adaptive())
    {
        this->evals_ = 1;
        return integrate_
        (
            dx*this->eqnPtr_->fx(0.5*(X0 + X1), li),
            X0,
            X1,
            this->tolerance_,
            li
        );
    }

    dx /= scalar(this->nIntervals_);
    scalar x12 = X0 + dx/2.0;
    Type res(dx*this->eqnPtr_->fx(x12, li));

    for (label i = 1; i < this->nIntervals_; i++)
    {
        x12 += dx;
        res = res + dx*this->eqnPtr_->fx(x12, li);
    }
    this->evals_ = this->nIntervals_;
    return res;
}
// ************************************************************************* //
