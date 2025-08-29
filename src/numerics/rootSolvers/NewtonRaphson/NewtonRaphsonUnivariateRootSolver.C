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

#include "NewtonRaphsonUnivariateRootSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace rootSolvers
{
namespace univariate
{
    defineTypeNameAndDebug(NewtonRaphson, 0);
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        NewtonRaphson,
        dictionaryOne
    );
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        NewtonRaphson,
        dictionaryTwo
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rootSolvers::univariate::NewtonRaphson::NewtonRaphson
(
    const scalarMultivariateEquation& eqn,
    const dictionary& dict
)
:
    univariateRootSolver(eqn, dict)
{}


Foam::rootSolvers::univariate::NewtonRaphson::NewtonRaphson
(
    const scalarMultivariateEquation& eqn,
    const NewtonRaphson& solver
)
:
    univariateRootSolver(eqn, solver)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rootSolvers::univariate::NewtonRaphson::~NewtonRaphson()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::rootSolvers::univariate::NewtonRaphson::findRoot
(
    const scalar x0,
    const scalar x1,
    const scalar x2,
    const label li
) const
{
    initialise(x0);
    scalar xOld = x0;
    scalar xNew = x0;
    scalar y = eqn_.fx(xOld, li);

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        xNew = xOld - y/stabilise(eqn_.dfdx(xOld, li), small);
        eqn_.limit(xNew);
        y = eqn_.fx(xNew, li);

        if (converged(xNew, xOld))
        {
            break;
        }

        xOld = xNew;
        printStepInformation(xNew);
    }

    return printFinalInformation(xNew);
}

// ************************************************************************* //
