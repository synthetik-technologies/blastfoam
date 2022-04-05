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

#include "multivariateEquation.H"

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::multivariateEquation<Type>::~multivariateEquation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
bool Foam::multivariateEquation<Type>::containsRoot
(
    const List<Type>& y0s,
    const List<Type>& y1s
) const
{
    forAll(y0s, i)
    {
        for (label cmpti = 0; cmpti < this->nVar(); cmpti++)
        {
            if ((component(y0s[i], cmpti)*component(y1s[i], cmpti)) > 0)
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
bool Foam::multivariateEquation<Type>::containsRoot(const label li) const
{
    List<Type> fxLow(nEqns());
    List<Type> fxHigh(nEqns());
    this->FX(lowerLimits(), li, fxLow);
    this->FX(upperLimits(), li, fxHigh);
    return containsRoot(fxLow, fxHigh);
}


// ************************************************************************* //
