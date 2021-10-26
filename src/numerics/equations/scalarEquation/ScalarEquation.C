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

Foam::ScalarEquation::ScalarEquation
(
    const label nVar,
    const scalarField& lowerLimit,
    const scalarField& upperLimit
)
:
    Equation<scalarField, scalar>(nVar, 1, lowerLimit, upperLimit)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ScalarEquation::~ScalarEquation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::ScalarEquation::nVar() const
{
    return Equation<scalarField, scalar>::nVar();
}


Foam::label Foam::ScalarEquation::nEqns() const
{
    return Equation<scalarField, scalar>::nEqns();
}


Foam::tmp<Foam::scalarField>
Foam::ScalarEquation::lowerLimits() const
{
    return Equation<scalarField, scalar>::lower();
}


Foam::tmp<Foam::scalarField>
Foam::ScalarEquation::upperLimits() const
{
    return Equation<scalarField, scalar>::upper();
}


Foam::tmp<Foam::scalarField>
Foam::ScalarEquation::dX() const
{
    return Equation<scalarField, scalar>::dx();
}


void Foam::ScalarEquation::setDX(const scalarField& newDx) const
{
    this->dx_ = newDx;
}


void Foam::ScalarEquation::limit(scalarField& x) const
{
    Equation<scalarField, scalar>::limit(x);
}

void Foam::ScalarEquation::calculateGradient
(
    const scalarField& x0,
    const label li,
    const scalar& fx0,
    scalarField& grad
) const
{
    scalar fx1;
    scalarField x1(x0);
    const scalarField dx(this->dX());
    const Equation<scalarField, scalar>& eqn(*this);
    for (label cmpti = 0; cmpti < this->nVar_; cmpti++)
    {
        x1 = x0;
        x1[cmpti] += dx[cmpti];
        eqn.f(x1, li, fx1);
        grad[cmpti] = (fx1 - fx0)/dx[cmpti];
    }
}


void Foam::ScalarEquation::gradient
(
    const scalarField& x0,
    const label li,
    scalar& fx0,
    scalarField& grad
) const
{
    const Equation<scalarField, scalar>& eqn(*this);
    eqn.f(x0, li, fx0);
    calculateGradient(x0, li, fx0, grad);
}


// ************************************************************************* //
