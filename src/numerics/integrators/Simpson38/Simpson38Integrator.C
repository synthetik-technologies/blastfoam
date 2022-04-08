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

#include "Simpson38Integrator.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Simpson38Integrator<Type>::Simpson38Integrator
(
    const equationType& eqn,
    const dictionary& dict
)
:
    Integrator<Type>(eqn, dict)
{}


template<class Type>
Foam::Simpson38Integrator<Type>::Simpson38Integrator
(
    const equationType& eqn,
    const integrator& inter
)
:
    Integrator<Type>(eqn, inter)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Simpson38Integrator<Type>::integrate_
(
    const scalar dx,
    const Type& f0,
    const Type& f13,
    const Type& f23,
    const Type& f1
) const
{
    return dx/8.0*(f0 + 3.0*(f13 + f23) + f1);
}


template<class Type>
Type Foam::Simpson38Integrator<Type>::integrate_
(
    const Type& Q,
    const scalar x0,
    const scalar x1,
    const Type& f0,
    const Type& f13,
    const Type& f23,
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

    this->intervals_++;
    const scalar x16 = x0 + 1.0/6.0*dx;
    const scalar x12 = x0 + 0.5*dx;
    const scalar x56 = x0 + 5.0/6.0*dx;

    const Type f16(this->eqnPtr_->fx(x16, li));
    const Type f12(this->eqnPtr_->fx(x12, li));
    const Type f56(this->eqnPtr_->fx(x56, li));
    this->evals_ += 3;

    const Type fx0(integrate_(x12 - x0, f0, f16, f13, f12));
    const Type fx1(integrate_(x1 - x12, f12, f23, f56, f1));
    const Type fx(fx0 + fx1);
    if (this->converged(fx, Q, dx, tol))
    {
        return fx;
    }
    else
    {
        return
            integrate_(fx0, x0, x12, f0, f16, f13, f12, tol/2.0, li)
          + integrate_(fx1, x12, x1, f12, f23, f56, f1, tol/2.0, li);
    }
}


template<class Type>
Type Foam::Simpson38Integrator<Type>::integrate
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
        const Type f13(this->eqnPtr_->fx(X0 + dx/3.0, li));
        const Type f23(this->eqnPtr_->fx(X0 + 2.0*dx/3.0, li));
        const Type f1(this->eqnPtr_->fx(X1, li));
        const Type Q(integrate_(dx, f0, f13, f23, f1));

        this->evals_ = 4;
        return integrate_
        (
            Q,
            X0, X1,
            f0, f13, f23, f1,
            this->tolerance_,
            li
        );
    }

    dx /= scalar(this->nIntervals_);
    scalar x13 = X0 + dx/3.0;
    scalar x23 = X0 + 2.0*dx/3.0;
    scalar x1 = X0 + dx;

    Type f0(this->eqnPtr_->fx(X0, li));
    Type f1(this->eqnPtr_->fx(x1, li));

    Type res
    (
        integrate_
        (
            dx,
            f0,
            this->eqnPtr_->fx(x13, li),
            this->eqnPtr_->fx(x23, li),
            f1
        )
    );
    for (label i = 1; i < this->nIntervals_; i++)
    {
        x13 += dx;
        x23 += dx;
        x1 += dx;

        f0 = f1;
        f1 = this->eqnPtr_->fx(x1, li);

        res = res
          + integrate_
            (
                dx,
                f0,
                this->eqnPtr_->fx(x13, li),
                this->eqnPtr_->fx(x23, li),
                f1
            );
    }
    this->evals_ = 3*this->nIntervals_ + 1;
    return res;
}

// ************************************************************************* //
