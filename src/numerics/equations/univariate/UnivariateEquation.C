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

#include "UnivariateEquation.H"
#include "adaptiveTypes.H"

// * * * * * * * * * * * * * * Static Data Functions * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::univariateEquation<Type>>
Foam::univariateEquation<Type>::New
(
    const dictionary& dict
)
{
    return New(dict.lookup<word>("type"), dict);
}

template<class Type>
Foam::autoPtr<Foam::univariateEquation<Type>>
Foam::univariateEquation<Type>::New
(
    const word& type,
    const dictionary& dict
)
{

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << typeName << " type "
            << type <<  nl << nl
            << "Valid " << typeName << " types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    return cstrIter()(dict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::UnivariateEquation<Type>::UnivariateEquation
(
    const label nVar,
    const string& eqnString
)
:
    lowerLimits_(nVar, -great),
    upperLimits_(nVar, great),
    nVar_(nVar),
    dX_(nVar_, 1e-6)
{}


template<class Type>
Foam::UnivariateEquation<Type>::UnivariateEquation
(
    const scalarList& lowerLimits,
    const scalarList& upperLimits,
    const string& eqnString
)
:
    univariateEquation<Type>(eqnString),
    lowerLimits_(lowerLimits),
    upperLimits_(upperLimits),
    nVar_(lowerLimits.size()),
    dX_(nVar_, 1e-6)
{}


template<class Type>
Foam::UnivariateEquation<Type>::UnivariateEquation
(
    const scalarList& lowerLimits,
    const scalarList& upperLimits,
    const dictionary& dict,
    const string& eqnString
)
:
    univariateEquation<Type>(eqnString, dict),
    lowerLimits_(dict.lookupOrDefault("lowerBounds", lowerLimits)),
    upperLimits_(dict.lookupOrDefault("upperBounds", upperLimits)),
    nVar_(lowerLimits.size()),
    dX_
    (
        dict.found("dx")
      ? scalarList(nVar_, dict.lookup<scalar>("dx"))
      : dict.lookupOrDefault<scalarList>("dX", scalarList(nVar_, 1e-6))
    )
{}


template<class Type>
Foam::UnivariateEquation<Type>::UnivariateEquation(const dictionary& dict)
:
    UnivariateEquation<Type>
    (
        dict.lookup("lowerBounds"),
        dict.lookup("upperBounds"),
        dict
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::UnivariateEquation<Type>::~UnivariateEquation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::UnivariateEquation<Type>::calculateGradient
(
    const typename univariateEquation<Type>::VarType& x0,
    const label li,
    List<Type>& grad
) const
{
    calculateGradient(this->fX(x0, li), x0, li, grad);
}


template<class Type>
void Foam::UnivariateEquation<Type>::calculateGradient
(
    const Type& fx0,
    const typename univariateEquation<Type>::VarType& x0,
    const label li,
    List<Type>& grad
) const
{
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
    const typename univariateEquation<Type>::VarType& x,
    const label li,
    List<Type>& fx
) const
{
    fx[0] = this->fX(x, li);
}


template<class Type>
void Foam::UnivariateEquation<Type>::dfdX
(
    const typename univariateEquation<Type>::VarType& x,
    const label li,
    List<Type>& dfdx
) const
{
    calculateGradient(x, li, dfdx);
}


template<class Type>
void Foam::UnivariateEquation<Type>::jacobian
(
    const typename univariateEquation<Type>::VarType& x,
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
    const UList<Type>& y0s,
    const UList<Type>& y1s
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
                << "Solution of component " << cmpti
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
