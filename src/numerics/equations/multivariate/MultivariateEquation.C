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

#include "MultivariateEquation.H"
#include "adaptiveTypes.H"

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
    lowerLimits_(lowerLimits),
    upperLimits_(upperLimits),
    nVar_(lowerLimits.size()),
    nEqns_(nEqns),
    dX_(nVar_, 1e-6)
{}


template<class Type>
Foam::MultivariateEquation<Type>::MultivariateEquation
(
    const List<string>& eqnStrings,
    const label nEqns,
    const scalarList& lowerLimits,
    const scalarList& upperLimits
)
:
    multivariateEquation<Type>(eqnStrings),
    lowerLimits_(lowerLimits),
    upperLimits_(upperLimits),
    nVar_(lowerLimits.size()),
    nEqns_(nEqns),
    dX_(nVar_, 1e-6)
{}


template<class Type>
Foam::MultivariateEquation<Type>::MultivariateEquation(const dictionary& dict)
:
    multivariateEquation<Type>(dict),
    lowerLimits_(dict.lookup("lowerBounds")),
    upperLimits_(dict.lookup("upperBounds")),
    nVar_(lowerLimits_.size()),
    nEqns_(dict.lookup<label>("nEquations")),
    dX_(dict.lookupOrDefault("dx", scalarList(nVar_, 1e-6)))
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


template<class Type>
bool Foam::MultivariateEquation<Type>::containsRoot
(
    const List<Type>& y0s,
    const List<Type>& y1s
) const
{
    forAll(y0s, i)
    {
        for (label cmpti = 0; cmpti < this->nVar(); cmpti++)
        {
            if
            (
                adaptiveError::cmpt<Type>(y0s[i], cmpti)
               *adaptiveError::cmpt<Type>(y1s[i], cmpti)
              > 0
            )
            {
                #ifdef FULLDEBUG
                FatalErrorInFunction
                    << "Solution of component " << i
                    << " is not bracked in "
                    << "(" << lowerLimits()
                    << ","<< upperLimits() << ")" << endl
                    << abort(FatalError);
                #endif
                return false;
            }
        }
    }
    return true;
}


template<class Type>
bool Foam::MultivariateEquation<Type>::containsRoot(const label li) const
{
    List<Type> fxLow(nEqns());
    List<Type> fxHigh(nEqns());
    this->FX(lowerLimits(), li, fxLow);
    this->FX(upperLimits(), li, fxHigh);
    return containsRoot(fxLow, fxHigh);
}


// ************************************************************************* //
