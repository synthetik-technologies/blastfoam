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

#include "integrator.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(integrator, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::integrator::integrator(const dictionary& dict)
:
    integratorBase(dict),
    tolerance_(dict.lookupOrDefault<scalar>("tolerance", 1e-6)),
    maxSplits_(dict.lookupOrDefault<label>("maxSplits", 10)),
    nIntervals_(dict.lookupOrDefault<label>("nIntervals", 10))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::integrator::reset(const scalar dx) const
{
    if (adaptive())
    {
        minDx_ = mag(dx)/scalar(pow(2, maxSplits_));
        intervals_ = 1;
    }
    else
    {
        intervals_ = nIntervals_;
    }
}

// ************************************************************************* //
