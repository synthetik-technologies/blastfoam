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

#include "TrapezoidalMultivariateIntegrator.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::TrapezoidalMultivariateIntegrator<Type>::TrapezoidalMultivariateIntegrator
(
    const equationType& eqn,
    const dictionary& dict
)
:
    MultivariateIntegrator<Type>(eqn, dict)
{}


template<class Type>
Foam::TrapezoidalMultivariateIntegrator<Type>::TrapezoidalMultivariateIntegrator
(
    const equationType& eqn,
    const multivariateIntegrator& inter
)
:
    MultivariateIntegrator<Type>(eqn, inter)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::TrapezoidalMultivariateIntegrator<Type>::integrateFunc
(
    const scalarField& x0,
    const scalarField& x1,
    const label li
) const
{
    scalarField x(x0);
    scalarField dx(x1 - x0);
    this->evals_++;
    Type fx(this->eqnPtr_->fX(x, li));

    addCorners(0, dx, li, x, fx);

    return 1.0/pow(2.0, x0.size())*fx;
}


template<class Type>
void Foam::TrapezoidalMultivariateIntegrator<Type>::addCorners
(
    const label diri,
    const scalarField& dx,
    const label li,
    scalarField& x,
    Type& fx
) const
{
    if (diri < dx.size())
    {
        addCorners(diri+1, dx, li, x, fx);

        const scalar xOrig(x[diri]);
        x[diri] += dx[diri];
        this->evals_++;
        fx = fx + this->eqnPtr_->fX(x, li);
        addCorners(diri+1, dx, li, x, fx);
        x[diri] = xOrig;
    }
}

template<class Type>
Type Foam::TrapezoidalMultivariateIntegrator<Type>::integrate
(
    const scalarField& x0,
    const scalarField& x1,
    const label li
) const
{
    scalarField dX(x1 - x0);
    this->reset(dX);

    scalarField xs(0.5*(x1 + x0));
    scalar dx = 1.0;
    forAll(x1, i)
    {
        dx *= x1[i] - x0[i];
    }
    Type Q(dx*integrateFunc(x0, x1, li));

    if (this->adaptive())
    {
        return this->integrate_(Q, 0, x0, x1, this->tolerance_, li);
    }
    return Q;
}


// ************************************************************************* //
