/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
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
namespace rootSolvers
{
namespace univariate
{
    defineTypeNameAndDebug(falsePoint, 0);
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        falsePoint,
        dictionaryZero
    );
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        falsePoint,
        dictionaryOne
    );
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        falsePoint,
        dictionaryTwo
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rootSolvers::univariate::falsePoint::falsePoint
(
    const scalarMultivariateEquation& eqn,
    const dictionary& dict
)
:
    univariateRootSolver(eqn, dict)
{}


Foam::rootSolvers::univariate::falsePoint::falsePoint
(
    const scalarMultivariateEquation& eqn,
    const falsePoint& solver
)
:
    univariateRootSolver(eqn, solver)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rootSolvers::univariate::falsePoint::~falsePoint()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::rootSolvers::univariate::falsePoint::findRoot
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
        xNew = (xHigh*yLow - xLow*yHigh)/stabilise(yLow - yHigh, small);
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
        if (converged(xHigh, xLow))
        {
            break;
        }


        printStepInformation(xNew);
    }

    return printFinalInformation(xNew);
}

// ************************************************************************* //
