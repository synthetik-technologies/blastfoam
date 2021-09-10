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

#include "multivariateRootSystem.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::multivariateRootSystem::checkLimits() const
{
    if (xMins_.size() != nEqns() || xMaxs_.size() != nEqns())
    {
        FatalErrorInFunction
            << "Limits have not been set, but are required for the " << nl
            << "requested root solver." << endl
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multivariateRootSystem::multivariateRootSystem()
{}


Foam::multivariateRootSystem::multivariateRootSystem
(
    const scalarField& xMins,
    const scalarField& xMaxs
)
:
    xMins_(xMins),
    xMaxs_(xMaxs)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multivariateRootSystem::~multivariateRootSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::multivariateRootSystem::checkConditions
(
    const scalarField& y0s,
    const scalarField& y1s
) const
{
    forAll(y0s, i)
    {
        if (y0s[i]*y1s[i] > 0)
        {
            FatalErrorInFunction
                << "Solution of component " << i << " is not bracked in "
                << "(" << xMins_[i] << ","<< xMaxs_[i] << ")" << endl
                << abort(FatalError);
        }
    }
}


void Foam::multivariateRootSystem::checkConditions(const label li) const
{
    scalarField y0s(nEqns());
    scalarField y1s(nEqns());
    f(xMins_, li, y0s);
    f(xMaxs_, li, y1s);
    checkConditions(y0s, y1s);
}

// ************************************************************************* //
