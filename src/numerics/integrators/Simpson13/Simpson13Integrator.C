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

#include "Simpson13Integrator.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Simpson13Integrator<Type>::Simpson13Integrator
(
    const equationType& eqn,
    const dictionary& dict
)
:
    Integrator<Type>(eqn, dict)
{}


template<class Type>
Foam::Simpson13Integrator<Type>::Simpson13Integrator
(
    const equationType& eqn,
    const integrator& inter
)
:
    Integrator<Type>(eqn, inter)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Simpson13Integrator<Type>::integrate_
(
    const scalar dx,
    const Type& f0,
    const Type& fm,
    const Type& f1
) const
{
    return dx/6.0*(f0 + 4.0*fm + f1);
}

template<class Type>
Type Foam::Simpson13Integrator<Type>::integrate_
(
    const Type& Q,
    const scalar x0,
    const scalar x1,
    const Type& f0,
    const Type& fm,
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

    const scalar xm = 0.5*(x0 + x1);
    const scalar x0m = 0.5*(x0 + xm);
    const scalar xm1 = 0.5*(xm + x1);

    const Type f0m(this->eqnPtr_->fx(x0m, li));
    const Type fm1(this->eqnPtr_->fx(xm1, li));
    this->evals_ += 2;

    const Type fx0(integrate_(xm - x0, f0, f0m, fm));
    const Type fx1(integrate_(x1 - xm, fm, fm1, f1));
    const Type fx(fx0 + fx1);
    if (this->converged(fx, Q, dx, tol))
    {
        return fx;
    }
    else
    {
        this->intervals_++;
        return
            integrate_(fx0, x0, xm, f0, f0m, fm, tol/2.0, li)
          + integrate_(fx1, xm, x1, fm, fm1, f1, tol/2.0, li);
    }
}


template<class Type>
Type Foam::Simpson13Integrator<Type>::integrate
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
        const Type f12(this->eqnPtr_->fx(X0 + 0.5*dx, li));
        const Type f1(this->eqnPtr_->fx(X1, li));

        this->evals_ = 3;
        return integrate_
        (
            integrate_(dx, f0, f12, f1),
            X0, X1,
            f0, f12, f1,
            this->tolerance_,
            li
        );
    }

    dx /= scalar(this->nIntervals_);
    scalar x12 = X0 + 0.5*dx;
    scalar x1 = X0 + dx;
    Type f0(this->eqnPtr_->fx(x1, li));
    Type f1(this->eqnPtr_->fx(x1, li));
    Type res
    (
        integrate_
        (
            dx,
            f0,
            this->eqnPtr_->fx(x12, li),
            f1
        )
    );
    for (label i = 1; i < this->nIntervals_; i++)
    {
        x12 += dx;
        x1 += dx;

        f0 = f1;
        f1 = this->eqnPtr_->fx(x1, li);

        res = res
          + integrate_
            (
                dx,
                f0,
                this->eqnPtr_->fx(x12, li),
                f1
            );
    }
    this->evals_ = 2*this->nIntervals_ + 1;
    return res;
}

// ************************************************************************* //
