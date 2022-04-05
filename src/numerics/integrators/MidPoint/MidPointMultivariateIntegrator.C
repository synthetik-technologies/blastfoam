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

#include "MidPointMultivariateIntegrator.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::MidPointMultivariateIntegrator<Type>::MidPointMultivariateIntegrator
(
    const equationType& eqn,
    const dictionary& dict
)
:
    MultivariateIntegrator<Type>(eqn, dict)
{}


template<class Type>
Foam::MidPointMultivariateIntegrator<Type>::MidPointMultivariateIntegrator
(
    const equationType& eqn,
    const multivariateIntegrator& inter
)
:
    MultivariateIntegrator<Type>(eqn, inter)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::MidPointMultivariateIntegrator<Type>::integrateFunc
(
    const scalarField& x0,
    const scalarField& x1,
    const label li
) const
{
    this->evals_++;
    scalarField xs(0.5*(x0 + x1));
    return this->eqnPtr_->fX(xs, li);
}


template<class Type>
Type Foam::MidPointMultivariateIntegrator<Type>::integrate
(
    const scalarField& x0,
    const scalarField& x1,
    const label li
) const
{
    scalarField dX(x1 - x0);
    this->reset(dX);

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
