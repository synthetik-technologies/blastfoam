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

#include "Equation.H"
#include "adaptiveTypes.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Equation<Type>::Equation
(
    const scalar lower,
    const scalar upper
)
:
    lower_(lower),
    upper_(upper)
{}


template<class Type>
Foam::Equation<Type>::Equation
(
    const string& name,
    const scalar lower,
    const scalar upper
)
:
    equation<Type>(name),
    lower_(lower),
    upper_(upper)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Equation<Type>::~Equation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Type Foam::Equation<Type>::fX
(
    const scalarList& x,
    const label li
) const
{
    return this->fx(x[0], li);
}


template<class Type>
void Foam::Equation<Type>::FX
(
    const scalarList& x,
    const label li,
    List<Type>& fx
) const
{
    fx[0] = this->fx(x[0], li);
}


template<class Type>
void Foam::Equation<Type>::dfdX
(
    const scalarList& x,
    const label li,
    List<Type>& dfdx
) const
{
    dfdx[0] = this->dfdx(x[0], li);
}


template<class Type>
void Foam::Equation<Type>::jacobian
(
    const scalarList& x,
    const label li,
    List<Type>& fx,
    RectangularMatrix<Type>& J
) const
{
    fx[0] = this->fx(x[0], li);
    J(0, 0) = this->dfdx(x[0], li);
}


template<class Type>
bool Foam::Equation<Type>::containsRoot(const label li) const
{
    return containsRoot(this->fx(lower_, li), this->fx(upper_, li));
}


template<class Type>
bool Foam::Equation<Type>::containsRoot
(
    const Type& y0,
    const Type& y1
) const
{
    for (label cmpti = 0; cmpti < adaptiveError::nCmpts<Type>(); cmpti++)
    {
        if
        (
            adaptiveError::cmpt<Type>(y0, cmpti)
           *adaptiveError::cmpt<Type>(y1, cmpti)
          > 0
        )
        {
            #ifdef FULLDEBUG
            FatalErrorInFunction
                << "Solution is not bracked:" << nl
                << "limits: (" << lower() << ","<< upper() << ")" << endl
                << "f(x0)=" << y0 << ", f(x1)=" << y1 << endl
                << abort(FatalError);
            #endif
            return false;
        }
    }
    return true;
}


template<class Type>
bool Foam::Equation<Type>::containsRoot
(
    const List<Type>& y0s,
    const List<Type>& y1s
) const
{
    return containsRoot(y0s[0], y1s[0]);
}


// ************************************************************************* //
