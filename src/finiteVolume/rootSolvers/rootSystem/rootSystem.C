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

#include "rootSystem.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rootSystem::rootSystem()
:
    xMin_(-great),
    xMax_(great)
{}


Foam::rootSystem::rootSystem(const scalar xMin, const scalar xMax)
:
    xMin_(xMin),
    xMax_(xMax)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rootSystem::~rootSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::rootSystem::checkConditions(const scalar y0, const scalar y1) const
{
    if (y0*y1 > 0)
    {
        FatalErrorInFunction
            << "Solution is not bracked in "
            << "(" << xMin_ << ","<< xMax_ << ")" << endl
            << abort(FatalError);
    }
}


void Foam::rootSystem::checkConditions(const label li) const
{
    checkConditions(f(xMin_, li), f(xMax_, li));
}


// ************************************************************************* //
