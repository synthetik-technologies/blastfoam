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

#include "falsePointUnivariateRootSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(falsePointUnivariateRootSolver, 0);
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        falsePointUnivariateRootSolver,
        dictionaryZero
    );
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        falsePointUnivariateRootSolver,
        dictionaryOne
    );
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        falsePointUnivariateRootSolver,
        dictionaryTwo
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::falsePointUnivariateRootSolver::falsePointUnivariateRootSolver
(
    const scalarMultivariateEquation& eqn,
    const dictionary& dict
)
:
    univariateRootSolver(eqn, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::falsePointUnivariateRootSolver::~falsePointUnivariateRootSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::falsePointUnivariateRootSolver::findRoot
(
    const scalar x0,
    const scalar x1,
    const scalar x2,
    const label li
) const
{
    initialise(x0);
    scalar xNew = x0;
    scalar xLow = x1;
    scalar xHigh = x2;
    scalar yLow = eqn_.fx(xLow, li);
    scalar yHigh = eqn_.fx(xHigh, li);

    if (!eqn_.containsRoot(yLow, yHigh))
    {
        return x0;
    }

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        xNew = (xHigh*yLow - xLow*yHigh)/stabilise(yLow - yHigh, 1e-10);
        eqn_.limit(xNew);
        scalar yNew = eqn_.fx(xNew, li);

        // Both on the same side so pick the smallest value
        if (yLow*yHigh > 0)
        {
            if (mag(yLow) < mag(yHigh))
            {
                xHigh = xNew;
                yHigh = yNew;
            }
            else
            {
                xLow = xNew;
                yLow = yNew;
            }
        }
        else if (yNew*yLow < 0)
        {
            xLow = xNew;
            yLow = yNew;
        }
        else
        {
            xHigh = xNew;
            yHigh = yNew;
        }
        if (converged(xHigh, xLow, yNew))
        {
            break;
        }


        printStepInformation(xNew);
    }

    return printFinalInformation(xNew);
}

// ************************************************************************* //
