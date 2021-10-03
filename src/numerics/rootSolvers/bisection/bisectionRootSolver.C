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
    addToRunTimeSelectionTable
    (
        rootSolver,
        bisectionRootSolver,
        dictionaryZero
    );
    addToRunTimeSelectionTable
    (
        rootSolver,
        bisectionRootSolver,
        dictionaryOne
    );
    addToRunTimeSelectionTable
    (
        rootSolver,
        bisectionRootSolver,
        dictionaryTwo
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bisectionRootSolver::bisectionRootSolver
(
    const scalarEquation& eqn,
    const dictionary& dict
)
:
    rootSolver(eqn, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::bisectionRootSolver::solve
(
    const scalar x0,
    const scalar x1,
    const scalar x2,
    const label li
) const
{
    scalar xMean = x0;
    scalar xLow = x1;
    scalar xHigh = x2;
    scalar y = eqn_.f(xMean, li);

    if (!eqn_.containsRoot(li))
    {
        return x0;
    }

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        if (y > 0)
        {
            xLow = xMean;
        }
        else
        {
            xHigh = xMean;
        }

        if (converged(y))
        {
            return xMean;
        }

        xMean = (xLow + xHigh)/2.0;
        y = eqn_.f(xMean, li);

    }
    WarningInFunction
        << "Could not converge to the given root." << endl;

    return xMean;
}

// ************************************************************************* //
