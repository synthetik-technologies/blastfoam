/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
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

#include "MultivariateIntegrator.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multivariateIntegrator, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multivariateIntegrator::multivariateIntegrator
(
    const label n,
    const dictionary& dict
)
:
    integratorBase(dict),
    tolerance_
    (
        dict.lookupOrDefault<List<scalar>>
        (
            "tolerance",
            List<scalar>(n, 1e-4)
        )
    ),
    absTolerance_
    (
        dict.lookupOrDefault<List<scalar>>
        (
            "absTolerance",
            List<scalar>(n, 1e-6)
        )
    ),
    maxSplits_
    (
        dict.lookupOrDefault<List<label>>
        (
            "maxSplits",
            List<label>(n, 5)
        )
    ),
    nIntervals_
    (
        dict.lookupOrDefault<List<label>>
        (
            "nIntervals",
            List<label>(n, 5)
        )
    ),
    intervals_(n, 0),
    minDx_(n, great)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multivariateIntegrator::reset(const List<scalar>& dx) const
{
    if (adaptive())
    {
        forAll(minDx_, i)
        {
            minDx_[i] = mag(dx[i])/scalar(pow(2, maxSplits_[i]));
            intervals_[i] = 1;
        }
    }
    else
    {
        intervals_ = nIntervals_;
    }
}

// ************************************************************************* //
