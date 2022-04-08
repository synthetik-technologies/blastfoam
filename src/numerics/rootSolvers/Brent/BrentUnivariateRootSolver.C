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

#include "BrentUnivariateRootSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BrentUnivariateRootSolver, 0);
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        BrentUnivariateRootSolver,
        dictionaryZero
    );
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        BrentUnivariateRootSolver,
        dictionaryOne
    );
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        BrentUnivariateRootSolver,
        dictionaryTwo
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BrentUnivariateRootSolver::BrentUnivariateRootSolver
(
    const scalarMultivariateEquation& eqn,
    const dictionary& dict
)
:
    univariateRootSolver(eqn, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::BrentUnivariateRootSolver::~BrentUnivariateRootSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::BrentUnivariateRootSolver::findRoot
(
    const scalar x,
    const scalar xLow,
    const scalar xHigh,
    const label li
) const
{
    initialise(x);
    scalar x0 = xLow;
    scalar x1 = xHigh;
    scalar xNew = x;
    scalar y0 = eqn_.fx(x0, li);
    scalar y1 = eqn_.fx(x1, li);

    if (!eqn_.containsRoot(y0, y1))
    {
        return x;
    }

    if (mag(y0) < mag(y1))
    {
        scalar xtmp = x0;
        scalar ytmp = y0;
        x0 = x1;
        y0 = y1;
        x1 = xtmp;
        y1 = ytmp;
    }

    scalar x2 = x0;
    scalar y2 = y0;
    scalar x3 = x2;
    bool bisection = true;

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        if (mag(y0 - y2)  < yTol() && mag(y1 - y2) < yTol())
        {
            xNew =
                x0*y1*y2/((y0 - y1)*(y0 - y2))
              + x1*y0*y2/((y1 - y0)*(y1 - y2))
              + x2*y0*y1/((y2 - y0)*(y2 - y1));
        }
        else
        {
            xNew = x1 - y1*(x1 - x0)/stabilise((y1 - y0), yTol());
        }
        eqn_.limit(xNew);

        // Use bisection method if satisfies the conditions.
        scalar delta = mag(xTol()*x1);
        scalar min1 = mag(xNew - x1);
        scalar min2 = mag(x1 - x2);
        scalar min3 = mag(x2 - x3);
        if
        (
            (xNew < (3.0*x0 + x1)/4.0 && xNew > x1)
         || (bisection && min1 >= min2/2.0)
         || (!bisection && min1 >= min3/2.0)
         || (bisection && min2 < delta)
         || (!bisection && min3 < delta)
        )
        {
            xNew = (x0 + x1)/2.0;
            bisection = true;
        }
        else
        {
            bisection = false;
        }

        scalar yNew = eqn_.fx(xNew, li);

        if (converged(x0, x1, yNew))
        {
            break;
        }

        x3 = x2;
        x2 = x1;

        if (y0*yNew < 0)
        {
            x1 = xNew;
            y1 = yNew;
        }
        else
        {
            x0 = xNew;
            y0 = yNew;
        }

        if (mag(y0) < mag(y1))
        {
            scalar xtmp = x0;
            scalar ytmp = y0;
            x0 = x1;
            y0 = y1;
            x1 = xtmp;
            y1 = ytmp;
        }
        printStepInformation(xNew);
    }
    return printFinalInformation(xNew);
}

// ************************************************************************* //
