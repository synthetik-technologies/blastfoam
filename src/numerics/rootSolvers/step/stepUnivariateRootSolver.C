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

#include "stepUnivariateRootSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace rootSolvers
{
namespace univariate
{
    defineTypeNameAndDebug(step, 0);
    addToRunTimeSelectionTable(univariateRootSolver, step, dictionaryZero);
    addToRunTimeSelectionTable(univariateRootSolver, step, dictionaryOne);
    addToRunTimeSelectionTable(univariateRootSolver, step, dictionaryTwo);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rootSolvers::univariate::step::step
(
    const scalarMultivariateEquation& eqn,
    const dictionary& dict
)
:
    univariateRootSolver(eqn, dict),
    dx_(dict.lookupOrDefault<scalar>("dx", -1.0)),
    f_(dict.lookupOrDefault<scalar>("f", 0.5))
{}


Foam::rootSolvers::univariate::step::step
(
    const scalarMultivariateEquation& eqn,
    const scalar dx
)
:
    univariateRootSolver(eqn, dictionary::null),
    dx_(dx),
    f_(0.5)
{}


Foam::rootSolvers::univariate::step::step
(
    const scalarMultivariateEquation& eqn,
    const step& solver
)
:
    univariateRootSolver(eqn, solver),
    dx_(solver.dx_),
    f_(solver.f_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rootSolvers::univariate::step::~step()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::rootSolvers::univariate::step::findRoot
(
    const scalar xm,
    const scalar x0,
    const scalar x1,
    const label li
) const
{
    initialise(xm);
    scalar x = x0;
    scalar dx = dx_ > 0 ? dx_ : (x1 - x0)/ceil(maxSteps_/100);
    scalar yLower = eqn_.fx(x0, li);

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        // New position
        scalar xNew = x + dx;

        // Limit position based on bounds
        if (xNew < x0)
        {
            xNew = x0;
        }
        else if (xNew > x1)
        {
            xNew = x1;
        }

        // Update dx (if limited)
        dx = xNew - x;

        // Value at new position
        const scalar y = eqn_.fx(xNew, li);

        // If root between x and xNew scale step size and retry
        if (y*yLower < 0)
        {
            dx *= f_;
        }

        // Advance to the next position
        else
        {
            x = xNew;
        }

        if (converged(dx))
        {
            break;
        }
        printStepInformation(x);

    }

    return printFinalInformation(x);
}

// ************************************************************************* //
