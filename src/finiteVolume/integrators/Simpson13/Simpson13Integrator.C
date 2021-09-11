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

#include "Simpson13Integrator.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Simpson13Integrator<Type>::Simpson13Integrator
(
    const Equation<Type>& eqn,
    const dictionary& dict
)
:
    Integrator<Type>(eqn, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Simpson13Integrator<Type>::integrate
(
    const scalar x0,
    const scalar x1,
    const label li
) const
{
    scalar dx = x1 - x0;
    if (dx < small)
    {
        return dx*this->eqn_.f(x0, li);
    }
    dx /= scalar(this->nSteps_);

    scalar a = x0;
    scalar b = a + dx;
    Type res
    (
        this->eqn_.f(a, li)
      + 4.0*this->eqn_.f(0.5*(a + b), li)
      + 2.0*this->eqn_.f(b, li)
    );
    for (label i = 1; i < this->nSteps_-1; i++)
    {
        a += dx;
        b += dx;
        res += 4.0*this->eqn_.f(0.5*(a + b), li) + 2.0*this->eqn_.f(b, li);
    }
    a += dx;
    b += dx;
    res += 4.0*this->eqn_.f(0.5*(a + b), li) + this->eqn_.f(b, li);
    return dx/6.0*res;
}

// ************************************************************************* //
