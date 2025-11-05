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

#include "bisectionUnivariateRootSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace rootSolvers
{
namespace univariate
{
    defineTypeNameAndDebug(bisection, 0);
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        bisection,
        dictionaryZero
    );
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        bisection,
        dictionaryOne
    );
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        bisection,
        dictionaryTwo
    );
}
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rootSolvers::univariate::bisection::bisection
(
    const scalarMultivariateEquation& eqn,
    const dictionary& dict
)
:
    univariateRootSolver(eqn, dict)
{}


Foam::rootSolvers::univariate::bisection::bisection
(
    const scalarMultivariateEquation& eqn,
    const bisection& solver
)
:
    univariateRootSolver(eqn, solver)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rootSolvers::univariate::bisection::~bisection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::rootSolvers::univariate::bisection::findRoot
(
    const scalar x0,
    const scalar x1,
    const scalar x2,
    const label li
) const
{
    initialise(x0);
    scalar xMean = x0;
    scalar xLow = x1;
    scalar xHigh = x2;
    scalar y = eqn_.fx(xMean, li);
    scalar yLow = eqn_.fx(xLow, li);
    scalar yHigh = eqn_.fx(xHigh, li);

    if (!eqn_.containsRoot(yLow, yHigh))
    {
        return x0;
    }

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        if (y*yHigh < 0)
        {
            xLow = xMean;
            yLow = y;
        }
        else
        {
            xHigh = xMean;
            yHigh = y;
        }

        if (converged(xLow, xHigh))
        {
            break;
        }

        xMean = (xLow + xHigh)/2.0;
        y = eqn_.fx(xMean, li);

        printStepInformation(xMean);

    }

    return printFinalInformation(xMean);
}

// ************************************************************************* //
