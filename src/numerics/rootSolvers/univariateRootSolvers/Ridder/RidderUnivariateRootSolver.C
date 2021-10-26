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

#include "RidderUnivariateRootSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RidderUnivariateRootSolver, 0);
    addToRunTimeSelectionTable
    (
        rootSolver,
        RidderUnivariateRootSolver,
        dictionaryUnivariate
    );
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        RidderUnivariateRootSolver,
        dictionaryZero
    );
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        RidderUnivariateRootSolver,
        dictionaryOne
    );
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        RidderUnivariateRootSolver,
        dictionaryTwo
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RidderUnivariateRootSolver::RidderUnivariateRootSolver
(
    const multivariateEquation<scalar>& eqn,
    const dictionary& dict
)
:
    univariateRootSolver(eqn, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RidderUnivariateRootSolver::~RidderUnivariateRootSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::RidderUnivariateRootSolver::findRoot
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
    scalar y0 = eqn_.fx(x0, li);
    scalar y1 = eqn_.fx(x1, li);

    if (!eqn_.containsRoot(y0, y1))
    {
        return x0;
    }

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        scalar xMean = 0.5*(x0 + x1);
        scalar yMean = eqn_.fx(xMean, li);

        xNew =
            xMean
          + (xMean - x0)*sign(y0 - y1)
           *yMean/sqrt(max(sqr(yMean) - y0*y1, small));

        if (converged(xNew - x0) || converged(xNew - x1))
        {
            break;
        }
        eqn_.limit(xNew);

        scalar yNew = eqn_.fx(xNew, li);
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
        printStepInformation(xNew);
    }

    return printFinalInformation(xNew);
}

// ************************************************************************* //
