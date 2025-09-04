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

#include "BrentUnivariateRootSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace rootSolvers
{
namespace univariate
{
    defineTypeNameAndDebug(Brent, 0);
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        Brent,
        dictionaryZero
    );
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        Brent,
        dictionaryOne
    );
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        Brent,
        dictionaryTwo
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rootSolvers::univariate::Brent::Brent
(
    const scalarMultivariateEquation& eqn,
    const dictionary& dict
)
:
    univariateRootSolver(eqn, dict)
{}


Foam::rootSolvers::univariate::Brent::Brent
(
    const scalarMultivariateEquation& eqn,
    const Brent& solver
)
:
    univariateRootSolver(eqn, solver)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rootSolvers::univariate::Brent::~Brent()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::rootSolvers::univariate::Brent::findRoot
(
    const scalar x,
    const scalar xLow,
    const scalar xHigh,
    const label li
) const
{
    initialise(x);
    scalar xa = xLow;
    scalar xb = xHigh;
    scalar ya = eqn_.fx(xa, li);
    scalar yb = eqn_.fx(xb, li);

    if (!eqn_.containsRoot(ya, yb))
    {
        return x;
    }

    if (mag(ya) < mag(yb))
    {
        Swap(xa, xb);
        Swap(ya, yb);
    }

    const scalar delta = relTolerance()*mag(x) + absTolerance();

    scalar xc = xa;
    scalar yc = ya;
    scalar xd = xc;
    scalar xNew = x;
    bool bisection = true;

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        if (mag(ya - yc) > small && mag(yb - yc) > small)
        {
            xNew =
                xa*yb*yc/((ya - yb)*(ya - yc))
              + xb*ya*yc/((yb - ya)*(yb - yc))
              + xc*ya*yb/((yc - ya)*(yc - yb));
        }
        else
        {
            xNew = xb - yb*(xb - xa)/stabilise(yb - ya, small);
        }

        eqn_.limit(xNew);

        // Use bisection method if satisfies the conditions.
        scalar xab = 0.25*(3.0*xa + xb);
        scalar min1 = mag(xNew - xb);
        scalar min2 = mag(xb - xc);
        scalar min3 = mag(xc - xd);

        if
        (
            ((xab - xNew)*(xNew - xb) < 0)
         || (bisection && min1 >= min2*0.5)
         || (!bisection && min1 >= min3*0.5)
         || (bisection && min2 < delta)
         || (!bisection && min3 < delta)
        )
        {
            xNew = (xa + xb)*0.5;
            bisection = true;
        }
        else
        {
            bisection = false;
        }

        if (converged(xa, xb))
        {
            break;
        }

        scalar yNew = eqn_.fx(xNew, li);


        xd = xc;
        xc = xb;
        yc = yb;

        if (ya*yNew < 0)
        {
            xb = xNew;
            yb = yNew;
        }
        else
        {
            xa = xNew;
            ya = yNew;
        }

        if (mag(ya) < mag(yb))
        {
            Swap(xa, xb);
            Swap(ya, yb);
        }
        printStepInformation(xNew);
    }
    return printFinalInformation(xNew);
}

// ************************************************************************* //
