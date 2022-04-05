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

#include "TrapezoidalIntegrator.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::TrapezoidalIntegrator<Type>::TrapezoidalIntegrator
(
    const equationType& eqn,
    const dictionary& dict
)
:
    Integrator<Type>(eqn, dict)
{}


template<class Type>
Foam::TrapezoidalIntegrator<Type>::TrapezoidalIntegrator
(
    const equationType& eqn,
    const integrator& inter
)
:
    Integrator<Type>(eqn, inter)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::TrapezoidalIntegrator<Type>::integrate_
(
    const scalar dx,
    const Type& f0,
    const Type& f1
) const
{
    return dx*0.5*(f0 + f1);
}


template<class Type>
Type Foam::TrapezoidalIntegrator<Type>::integrate_
(
    const Type& Q,
    const scalar x0,
    const scalar x1,
    const Type& f0,
    const Type& f1,
    const scalar tol,
    const label li
) const
{
    const scalar dx(x1 - x0);
    if (mag(dx) <= this->minDx_)
    {
        return Q;
    }

    const scalar x12 = x0 + 0.5*dx;
    const Type f12(this->eqnPtr_->fx(x12, li));
    this->evals_++;

    const Type fx0(integrate_(x12 - x0, f0, f12));
    const Type fx1(integrate_(x1 - x12, f12, f1));
    const Type fx(fx0 + fx1);
    if (this->converged(fx, Q, dx, tol))
    {
        return fx;
    }
    else
    {
        this->intervals_++;
        return
            integrate_(fx0, x0, x12, f0, f12, tol/2.0, li)
          + integrate_(fx1, x12, x1, f12, f1, tol/2.0, li);
    }
}


template<class Type>
Type Foam::TrapezoidalIntegrator<Type>::integrate
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
        const Type f0(this->eqnPtr_->fx(X0, li));
        const Type f1(this->eqnPtr_->fx(X1, li));

        this->evals_ = 2;
        return integrate_
        (
            integrate_(dx, f0, f1),
            X0, X1,
            f0, f1,
            this->tolerance_,
            li
        );
    }

    dx /= scalar(this->nIntervals_);
    scalar x1 = X0 + dx;
    Type f0(this->eqnPtr_->fx(x1, li));
    Type f1(this->eqnPtr_->fx(x1, li));
    Type res(integrate_(dx,f0, f1));
    for (label i = 1; i < this->nIntervals_; i++)
    {
        x1 += dx;
        f0 = f1;
        f1 = this->eqnPtr_->fx(x1, li);

        res = res + integrate_(dx, f0, f1);
    }
    this->evals_ = this->nIntervals_ + 1;
    return res;
}

// ************************************************************************* //
