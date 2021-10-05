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

#include "scalarEquation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scalarEquation::scalarEquation()
{}


Foam::scalarEquation::scalarEquation
(
    const scalar lowerLimit,
    const scalar upperLimit
)
:
    Equation<scalar>(lowerLimit, upperLimit)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::scalarEquation::~scalarEquation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::scalarEquation::containsRoot(const scalar y0, const scalar y1) const
{
    if (y0*y1 > 0)
    {
        #ifdef FULLDEBUG
        FatalErrorInFunction
            << "Solution is not bracked:" << nl
            << "limits: (" << lowerLimit_ << ","<< upperLimit_ << ")" << endl
            << "f(x0)=" << y0 << ", f(x1)=" << y1 << endl
            << abort(FatalError);
        #endif
        return false;
    }
    return true;
}


bool Foam::scalarEquation::containsRoot(const label li) const
{
    return containsRoot(f(lowerLimit_, li), f(upperLimit_, li));
}


// ************************************************************************* //
