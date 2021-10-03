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

#include "BrentRootSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BrentRootSolver, 0);
    addToRunTimeSelectionTable
    (
        rootSolver,
        BrentRootSolver,
        dictionaryZero
    );
    addToRunTimeSelectionTable
    (
        rootSolver,
        BrentRootSolver,
        dictionaryOne
    );
    addToRunTimeSelectionTable
    (
        rootSolver,
        BrentRootSolver,
        dictionaryTwo
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BrentRootSolver::BrentRootSolver
(
    const scalarEquation& eqn,
    const dictionary& dict
)
:
    rootSolver(eqn, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::BrentRootSolver::solve
(
    const scalar x,
    const label li
) const
{
    scalar x0 = eqn_.lower();
    scalar x1 = eqn_.upper();
    scalar xNew = x;
    scalar y0 = eqn_.f(x0, li);
    scalar y1 = eqn_.f(x1, li);

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
        if (converged(x1 - x0))
        {
            return x1;
        }

        if (mag(y0 - y2) > tolerance_ && mag(y1 - y2) > tolerance_)
        {
            xNew =
                x0*y1*y2/((y0 - y1)*(y0 - y2))
              + x1*y0*y2/((y1 - y0)*(y1 - y2))
              + x2*y0*y1/((y2 - y0)*(y2 - y1));
        }
        else
        {
            xNew = x1 - y1*(x1 - x0)/stabilise((y1 - y0), tolerance_);
        }

        // Use bisection method if satisfies the conditions.
        scalar delta = mag(tolerance_*x1);
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

        scalar yNew = eqn_.f(xNew, li);

        if (converged(yNew))
        {
            return xNew;
        }

        x3 = x2;
        x2 = x1;

        if (y0*y1 < 0)
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
    }
    WarningInFunction
        << "Could not converge to the given root." << endl;

    return xNew;
}

// ************************************************************************* //
