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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::MultivariateEquation<Type>::checkLimits() const
{
    if (lowerLimits_.size() != nEqns() || upperLimits_.size() != nEqns())
    {
        FatalErrorInFunction
            << "Limits have not been set, but are required for the " << nl
            << "requested root solver." << endl
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::MultivariateEquation<Type>::MultivariateEquation()
{}


template<class Type>
Foam::MultivariateEquation<Type>::MultivariateEquation
(
    const scalarField& lowerLimits,
    const scalarField& upperLimits
)
:
    lowerLimits_(lowerLimits),
    upperLimits_(upperLimits)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::MultivariateEquation<Type>::~MultivariateEquation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
bool Foam::MultivariateEquation<Type>::checkBounds(const scalarField& xs) const
{
    checkLimits();
    forAll(xs, i)
    {
        if (xs[i] < lowerLimits_[i] || xs[i] > upperLimits_[i])
        {
            #ifdef FULLDEBUG
            FatalErrorInFunction
                << "Request function evaluation is out of bounds." << nl
                << "lowerLimits: " << lowerLimits_ << endl
                << "upperLimits: " << upperLimits_ << endl
                << "x: " << xs << endl
                << abort(FatalError);
            #endif
            return false;
        }
    }
    return true;
}

// ************************************************************************* //
