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

#include "Simpson38Integrator.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Simpson38Integrator<Type>::Simpson38Integrator
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
Type Foam::Simpson38Integrator<Type>::integrate
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
    if (this->nSteps_ <= 1)
    {
        return
            dx/8.0
           *(
                this->eqnPtr_->fx(x0, li)
              + 3.0*this->eqnPtr_->fx((2.0*x0 + x1)/3.0, li)
              + 3.0*this->eqnPtr_->fx((x0 + 2.0*x1)/3.0, li)
              + this->eqnPtr_->fx(x1, li)
            );
    }
    dx /= scalar(this->nSteps_);

    Type res(this->eqnPtr_->fx(x0, li) + this->eqnPtr_->fx(x1, li));
    for (label i = 1; i <= this->nIntervals_; i++)
    {
        res +=
            3.0
           *(
                this->eqnPtr_->fx(x0 + dx*(3*i-2), li)
              + this->eqnPtr_->fx(x0 + dx*(3*i-1), li)
            );
    }
    for (label i = 1; i < this->nIntervals_; i++)
    {
        res += 2.0*this->eqnPtr_->fx(x0 + i*dx*3.0, li);
    }

    return dx*3.0/8.0*res;
}

// ************************************************************************* //
