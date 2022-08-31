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

#include "BooleIntegrator.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::BooleIntegrator<Type>::BooleIntegrator
(
    const equationType& eqn,
    const dictionary& dict
)
:
    Integrator<Type>(eqn, dict)
{}


template<class Type>
Foam::BooleIntegrator<Type>::BooleIntegrator
(
    const equationType& eqn,
    const integrator& inter
)
:
    Integrator<Type>(eqn, inter)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::BooleIntegrator<Type>::integrate_
(
    const scalar dx,
    const Type& f0,
    const Type& f14,
    const Type& f12,
    const Type& f34,
    const Type& f1
) const
{
    return dx/90.0*(7.0*(f0 + f1) + 32.0*(f14 + f34) + 12.0*f12);
}


template<class Type>
Type Foam::BooleIntegrator<Type>::integrate_
(
    const Type& Q,
    const scalar x0,
    const scalar x1,
    const Type& f0,
    const Type& f14,
    const Type& f12,
    const Type& f34,
    const Type& f1,
    const scalar tol,
    const label li
) const
{
    const scalar dx(x1 - x0);
    if (mag(dx) < this->minDx_)
    {
        return Q;
    }

    this->intervals_++;
    const scalar x18 = x0 + 0.125*dx;
    const scalar x38 = x0 + 0.375*dx;
    const scalar x12 = x0 + 0.5*dx;
    const scalar x58 = x0 + 0.625*dx;
    const scalar x78 = x0 + 0.875*dx;

    const Type f18(this->eqnPtr_->fx(x18, li));
    const Type f38(this->eqnPtr_->fx(x38, li));
    const Type f58(this->eqnPtr_->fx(x58, li));
    const Type f78(this->eqnPtr_->fx(x78, li));
    this->evals_ += 4;

    const Type fx0(integrate_(x12 - x0, f0, f18, f14, f38, f12));
    const Type fx1(integrate_(x1 - x12, f12, f58, f34, f78, f1));
    const Type fx(fx0 + fx1);
    if (this->converged(fx, Q, dx, tol))
    {
        return fx;
    }
    else
    {
        return
            integrate_(fx0, x0, x12, f0, f18, f14, f38, f12, tol/2.0, li)
          + integrate_(fx1, x12, x1, f12, f58, f34, f78, f1, tol/2.0, li);
    }
}


template<class Type>
Type Foam::BooleIntegrator<Type>::integrate
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
        const Type f14(this->eqnPtr_->fx(X0 + 0.25*dx, li));
        const Type f12(this->eqnPtr_->fx(X0 + 0.5*dx, li));
        const Type f34(this->eqnPtr_->fx(X0 + 0.75*dx, li));
        const Type f1(this->eqnPtr_->fx(X1, li));
        const Type Q(integrate_(dx, f0, f14, f12, f34, f1));

        this->evals_ = 5;
        return integrate_
        (
            Q,
            X0, X1,
            f0, f14, f12, f34, f1,
            this->tolerance_,
            li
        );
    }

    dx /= scalar(this->nIntervals_);
    scalar x14 = X0 + 0.25*dx;
    scalar x12 = X0 + 0.5*dx;
    scalar x34 = X0 + 0.75*dx;
    scalar x1 = X0 + dx;

    Type f0(this->eqnPtr_->fx(X0, li));
    Type f1(this->eqnPtr_->fx(x1, li));

    Type res
    (
        integrate_
        (
            dx,
            f0,
            this->eqnPtr_->fx(x14, li),
            this->eqnPtr_->fx(x12, li),
            this->eqnPtr_->fx(x34, li),
            f1
        )
    );
    for (label i = 1; i < this->nIntervals_; i++)
    {
        x14 += dx;
        x12 += dx;
        x34 += dx;
        x1 += dx;

        f0 = f1;
        f1 = this->eqnPtr_->fx(x1, li);

        res = res
          + integrate_
            (
                dx,
                f0,
                this->eqnPtr_->fx(x14, li),
                this->eqnPtr_->fx(x12, li),
                this->eqnPtr_->fx(x34, li),
                f1
            );
    }
    this->evals_ = 4*this->nIntervals_ + 1;
    return res;
}

// ************************************************************************* //
