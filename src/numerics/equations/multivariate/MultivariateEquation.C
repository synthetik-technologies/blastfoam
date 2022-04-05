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
    const scalarList& x0,
    const label li,
    const List<Type>& f0,
    RectangularMatrix<Type>& J
) const
{
    scalarList f1(nVar_);
    J.setSize(nEqns_, nVar_);
    for (label cmptj = 0; cmptj < nVar_; cmptj++)
    {
        scalarList x1(x0);
        x1[cmptj] += dX_[cmptj];
        this->FX(x1, li, f1);
        for (label cmpti = 0; cmpti < nEqns_; cmpti++)
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
    const scalarList& lowerLimits,
    const scalarList& upperLimits
)
:
    nVar_(lowerLimits.size()),
    nEqns_(nEqns),
    lowerLimits_(lowerLimits),
    upperLimits_(upperLimits),
    dX_(nVar_, 1e-6)
{}


template<class Type>
Foam::MultivariateEquation<Type>::MultivariateEquation
(
    const string& name,
    const label nEqns,
    const scalarList& lowerLimits,
    const scalarList& upperLimits
)
:
    multivariateEquation<Type>(name),
    nVar_(lowerLimits.size()),
    nEqns_(nEqns),
    lowerLimits_(lowerLimits),
    upperLimits_(upperLimits),
    dX_(nVar_, 1e-6)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::MultivariateEquation<Type>::~MultivariateEquation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::MultivariateEquation<Type>::jacobian
(
    const scalarList& x,
    const label li,
    List<Type>& fx,
    RectangularMatrix<Type>& J
) const
{
    fx.resize(nEqns_);
    this->FX(x, li, fx);
    calculateJacobian(x, li, fx, J);
}


// ************************************************************************* //
