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
    scalarField f1(this->nVar_);
    J.setSize(this->nEqns_, this->nVar_);
    const Equation<scalarField, Field<Type>>& eqns(*this);
    for (label cmptj = 0; cmptj < this->nVar_; cmptj++)
    {
        scalarField x1(x0);
        x1[cmptj] += this->dx(cmptj);
        eqns.f(x1, li, f1);
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

template<class Type>
Foam::label Foam::MultivariateEquation<Type>::nVar() const
{
    return Equation<scalarField, Field<Type>>::nVar();
}


template<class Type>
Foam::label Foam::MultivariateEquation<Type>::nEqns() const
{
    return Equation<scalarField, Field<Type>>::nEqns();
}


template<class Type>
Foam::tmp<Foam::scalarField>
Foam::MultivariateEquation<Type>::lowerLimits() const
{
    return Equation<scalarField, Field<Type>>::lower();
}


template<class Type>
Foam::tmp<Foam::scalarField>
Foam::MultivariateEquation<Type>::upperLimits() const
{
    return Equation<scalarField, Field<Type>>::lower();
}


template<class Type>
Foam::tmp<Foam::scalarField>
Foam::MultivariateEquation<Type>::dX() const
{
    return Equation<scalarField, Field<Type>>::dx();
}


template<class Type>
void Foam::MultivariateEquation<Type>::setDX
(
    const scalarField& newDx
) const
{
    this->dx_ = newDx;
}


template<class Type>
void Foam::MultivariateEquation<Type>::limit(scalarField& x) const
{
    Equation<scalarField, Field<Type>>::limit(x);
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
    const Equation<scalarField, Field<Type>>& eqns(*this);
    eqns.f(x, li, fx);
    calculateJacobian(x, li, fx, J);
}


// ************************************************************************* //
