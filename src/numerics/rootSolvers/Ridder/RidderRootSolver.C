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

#include "RidderRootSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RidderRootSolver, 0);
    addToRunTimeSelectionTable
    (
        rootSolver,
        RidderRootSolver,
        dictionaryZero
    );
    addToRunTimeSelectionTable
    (
        rootSolver,
        RidderRootSolver,
        dictionaryOne
    );
    addToRunTimeSelectionTable
    (
        rootSolver,
        RidderRootSolver,
        dictionaryTwo
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RidderRootSolver::RidderRootSolver
(
    const scalarEquation& eqn,
    const dictionary& dict
)
:
    rootSolver(eqn, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::RidderRootSolver::findRoot
(
    const scalar x,
    const scalar xLow,
    const scalar xHigh,
    const label li
) const
{
    scalar x0 = xLow;
    scalar x1 = xHigh;
    scalar xNew = x;
    scalar y0 = eqn_.f(x0, li);
    scalar y1 = eqn_.f(x1, li);

    if (!eqn_.containsRoot(y0, y1))
    {
        return x0;
    }

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        scalar xMean = 0.5*(x0 + x1);
        scalar yMean = eqn_.f(xMean, li);

        xNew =
            xMean
          + (xMean - x0)*sign(y0 - y1)
           *yMean/sqrt(max(sqr(yMean) - y0*y1, small));

        if (converged(xNew - x0) || converged(xNew - x1))
        {
            return xNew;
        }
        limit(xNew);

        scalar yNew = eqn_.f(xNew, li);
        if (converged(yNew))
        {
            return xNew;
        }

        if (yMean*yNew < 0)
        {
            x0 = xMean;
            y0 = yMean;
            x1 = xNew;
            y1 = yNew;
        }
        else if (y1*yNew < 0)
        {
            x0 = xNew;
            y0 = yNew;
        }
        else
        {
            x1 = xNew;
            y1 = yNew;
        }
    }
    printNoConvergence();

    return xNew;
}

// ************************************************************************* //
