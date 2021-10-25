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
    const equationType& eqn,
    const dictionary& dict
)
:
    Integrator<Type>(eqn, dict)
{
    this->setNIntervals(this->nIntervals_);
}


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
    if (mag(dx) < small)
    {
        return dx*this->eqnPtr_->fx(x0, li);
    }
    if (this->nSteps_ <= 2)
    {
        return
            dx/6.0
           *(
                this->eqnPtr_->fx(x0, li)
              + 4.0*this->eqnPtr_->fx(0.5*(x0 + x1), li)
              + this->eqnPtr_->fx(x1, li)
            );
    }
    label n = this->nSteps_/2.0;
    dx /= scalar(this->nSteps_);

    Type res(this->eqnPtr_->fx(x0, li) + this->eqnPtr_->fx(x1, li));
    for (label i = 1; i <= n; i++)
    {
        res += 4.0*this->eqnPtr_->fx(x0 + dx*(2*i - 1), li);
    }
    for (label i = 1; i < n; i++)
    {
        res += 2.0*this->eqnPtr_->fx(x0 + dx*2*i, li);
    }
    return dx/3.0*res;
}

// ************************************************************************* //
