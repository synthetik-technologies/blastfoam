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

#include "UnivariateEquation.H"
#include "adaptiveTypes.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::UnivariateEquation<Type>::UnivariateEquation
(
    const scalarList& lowerLimits,
    const scalarList& upperLimits
)
:
    nVar_(lowerLimits.size()),
    lowerLimits_(lowerLimits),
    upperLimits_(upperLimits),
    dX_(nVar_, 1e-6)
{}


template<class Type>
Foam::UnivariateEquation<Type>::UnivariateEquation
(
    const string& name,
    const scalarList& lowerLimits,
    const scalarList& upperLimits
)
:
    univariateEquation<Type>(name),
    nVar_(lowerLimits.size()),
    lowerLimits_(lowerLimits),
    upperLimits_(upperLimits),
    dX_(nVar_, 1e-6)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::UnivariateEquation<Type>::~UnivariateEquation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::UnivariateEquation<Type>::calculateGradient
(
    const scalarList& x0,
    const label li,
    List<Type>& grad
) const
{
    const Type fx0(this->fX(x0, li));
    scalarList x1(x0);
    for (label cmpti = 0; cmpti < nVar_; cmpti++)
    {
        x1[cmpti] += dX_[cmpti];
        grad[cmpti] = (this->fX(x1, li) - fx0)/dX_[cmpti];
        x1[cmpti] = x0[cmpti];
    }
}


template<class Type>
void Foam::UnivariateEquation<Type>::FX
(
    const scalarList& x,
    const label li,
    List<Type>& fx
) const
{
    fx[0] = this->fX(x, li);
}


template<class Type>
void Foam::UnivariateEquation<Type>::dfdX
(
    const scalarList& x,
    const label li,
    List<Type>& dfdx
) const
{
    calculateGradient(x, li, dfdx);
}


template<class Type>
void Foam::UnivariateEquation<Type>::jacobian
(
    const scalarList& x,
    const label li,
    List<Type>& fx,
    RectangularMatrix<Type>& J
) const
{
    fx.setSize(1);
    fx[0] = this->fX(x, li);

    J.setSize(1, x.size());
    List<Type> dfdx(x.size());
    this->dfdX(x, li, dfdx);
    forAll(dfdx, i)
    {
        J(0, i) = dfdx[i];
    }
}


template<class Type>
bool Foam::UnivariateEquation<Type>::containsRoot
(
    const List<Type>& y0s,
    const List<Type>& y1s
) const
{
    for (label cmpti = 0; cmpti < this->nVar(); cmpti++)
    {
        if
        (
            adaptiveError::cmpt<Type>(y0s[0], cmpti)
            *adaptiveError::cmpt<Type>(y1s[0], cmpti)
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
    return true;
}


template<class Type>
bool Foam::UnivariateEquation<Type>::containsRoot(const label li) const
{
    List<Type> fxLow(nEqns());
    List<Type> fxHigh(nEqns());
    this->FX(lowerLimits(), li, fxLow);
    this->FX(upperLimits(), li, fxHigh);
    return containsRoot(fxLow, fxHigh);
}


// ************************************************************************* //
