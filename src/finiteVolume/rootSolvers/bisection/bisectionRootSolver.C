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

#include "bisectionRootSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bisectionRootSolver, 0);
    addToRunTimeSelectionTable(rootSolver, bisectionRootSolver, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bisectionRootSolver::bisectionRootSolver
(
    const rootSystem& eqn,
    const dictionary& dict
)
:
    rootSolver(eqn, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::bisectionRootSolver::solve
(
    const scalar x0,
    const label li
) const
{
    scalar xMean = x0;
    scalar xLow = eqn_.lower();
    scalar xHigh = eqn_.upper();

    eqn_.checkConditions(li);

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        xMean = (xLow + xHigh)/2.0;
        scalar y = eqn_.f(xMean, li);

        if (y > 0)
        {
            xLow = xMean;
        }
        else
        {
            xHigh = xMean;
        }

        if (converged(xHigh - xLow) || converged(y))
        {
            return xMean;
        }

    }
    WarningInFunction
        << "Could not converge to the given root." << endl;

    return xMean;
}

// ************************************************************************* //
