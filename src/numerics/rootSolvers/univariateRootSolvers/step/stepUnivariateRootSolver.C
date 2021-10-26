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

#include "stepUnivariateRootSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(stepUnivariateRootSolver, 0);
    addToRunTimeSelectionTable
    (
        rootSolver,
        stepUnivariateRootSolver,
        dictionaryUnivariate
    );
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        stepUnivariateRootSolver,
        dictionaryZero
    );
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        stepUnivariateRootSolver,
        dictionaryOne
    );
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        stepUnivariateRootSolver,
        dictionaryTwo
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::stepUnivariateRootSolver::stepUnivariateRootSolver
(
    const multivariateEquation<scalar>& eqn,
    const dictionary& dict
)
:
    univariateRootSolver(eqn, dict),
    dx_
    (
        dict.lookupOrDefault<scalar>
        (
            "dx",
            (eqn_.upper() - eqn_.lower())
           /ceil(maxSteps_/10)
        )
    )
{}


Foam::stepUnivariateRootSolver::stepUnivariateRootSolver
(
    const equation& eqn,
    const scalar dx
)
:
    univariateRootSolver(eqn, dictionary()),
    dx_(dx)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::stepUnivariateRootSolver::~stepUnivariateRootSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::stepUnivariateRootSolver::findRoot
(
    const scalar x0,
    const scalar x1,
    const scalar x2,
    const label li
) const
{
    scalar x = x1;
    scalar dx = dx_;
    scalar yLower = eqn_.fx(x0, li);

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        scalar y = eqn_.fx(x + dx, li);
        if (y*yLower < 0)
        {
            dx /= 2.0;
        }
        else
        {
            x += dx;
        }
        if (converged(dx) || converged(y))
        {
            break;
        }
        printStepInformation(x);

    }

    return printFinalInformation(x);
}

// ************************************************************************* //
