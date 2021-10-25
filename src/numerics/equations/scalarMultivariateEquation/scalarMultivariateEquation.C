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

#include "scalarMultivariateEquation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scalarMultivariateEquation::scalarMultivariateEquation
(
    const label nEqns,
    const scalarField& lowerLimits,
    const scalarField& upperLimits
)
:
    MultivariateEquation<scalar>
    (
        nEqns,
        lowerLimits, 
        upperLimits
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::scalarMultivariateEquation::~scalarMultivariateEquation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::scalarMultivariateEquation::containsRoot
(
    const scalarField& y0s,
    const scalarField& y1s
) const
{
    forAll(y0s, i)
    {
        if (y0s[i]*y1s[i] > 0)
        {
            #ifdef FULLDEBUG
            FatalErrorInFunction
                << "Solution of component " << i << " is not bracked in "
                << "(" << lowerLimits_[i] << ","<< upperLimits_[i] << ")" << endl
                << abort(FatalError);
            #endif
            return false;
        }
    }
    return true;
}


bool Foam::scalarMultivariateEquation::containsRoot(const label li) const
{
    scalarField fxLow(nEqns());
    scalarField fxHigh(nEqns());
    this->f(lowerLimits_, li, fxLow);
    this->f(upperLimits_, li, fxHigh);
    return containsRoot(fxLow, fxHigh);
}

// ************************************************************************* //
