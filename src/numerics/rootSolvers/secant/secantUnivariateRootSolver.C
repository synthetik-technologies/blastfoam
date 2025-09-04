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

#include "secantUnivariateRootSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace rootSolvers
{
namespace univariate
{
    defineTypeNameAndDebug(secant, 0);
    addToRunTimeSelectionTable(univariateRootSolver, secant, dictionaryZero);
    addToRunTimeSelectionTable(univariateRootSolver, secant, dictionaryOne);
    addToRunTimeSelectionTable(univariateRootSolver, secant, dictionaryTwo);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rootSolvers::univariate::secant::secant
(
    const scalarMultivariateEquation& eqn,
    const dictionary& dict
)
:
    univariateRootSolver(eqn, dict)
{}


Foam::rootSolvers::univariate::secant::secant
(
    const scalarMultivariateEquation& eqn,
    const secant& solver
)
:
    univariateRootSolver(eqn, solver)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rootSolvers::univariate::secant::~secant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::rootSolvers::univariate::secant::findRoot
(
    const scalar x0,
    const scalar x1,
    const scalar x2,
    const label li
) const
{
    initialise(x0);
    scalar dx = relTolerance()*mag(x0) + absTolerance();
    scalar xPrev = x0 - dx;
    if (xPrev < eqn_.lower())
    {
        xPrev = x0 + dx;
    }

    scalar x = x0;
    scalar yPrev = eqn_.fx(xPrev, li);

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        scalar y = eqn_.fx(x, li);
        scalar xNew = x - y*(x - xPrev)/stabilise(y - yPrev, small);
        eqn_.limitChange(x, xNew, boundsFac_);

        xPrev = x;
        x = xNew;

        yPrev = y;

        if (converged(x, xPrev))
        {
            break;
        }
        printStepInformation(x);

    }

    return printFinalInformation(x);
}

// ************************************************************************* //
