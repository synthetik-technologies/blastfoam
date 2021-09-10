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

#include "falsePointRootSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(falsePointRootSolver, 0);
    addToRunTimeSelectionTable(rootSolver, falsePointRootSolver, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::falsePointRootSolver::falsePointRootSolver
(
    const rootSystem& eqn,
    const dictionary& dict
)
:
    rootSolver(eqn, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::falsePointRootSolver::solve
(
    const scalar x0,
    const label li
) const
{
    scalar xNew = x0;
    scalar xLow = eqn_.lower();
    scalar xHigh = eqn_.upper();
    scalar yLow = eqn_.f(xLow, li);
    scalar yHigh = eqn_.f(xHigh, li);

    eqn_.checkConditions(yLow, yHigh);

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        xNew = (xHigh*yLow - xLow*yHigh)/stabilise(yLow - yHigh, tolerance_);
        scalar yNew = eqn_.f(xNew, li);

        if (converged(yNew))
        {
            return xNew;
        }

        if (yNew > 0)
        {
            xLow = xNew;
            yLow = yNew;
        }
        else
        {
            xHigh = xNew;
            yHigh = yNew;
        }
        if (converged(xHigh - xLow))
        {
            return xNew;
        }
    }
    WarningInFunction
        << "Could not converge to the given root." << endl;

    return xNew;
}

// ************************************************************************* //
