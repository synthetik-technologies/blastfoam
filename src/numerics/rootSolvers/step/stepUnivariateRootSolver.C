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
    const scalarMultivariateEquation& eqn,
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
    ),
    f_(dict.lookupOrDefault<scalar>("f", 0.5))
{}


Foam::stepUnivariateRootSolver::stepUnivariateRootSolver
(
    const scalarMultivariateEquation& eqn,
    const scalar dx
)
:
    univariateRootSolver(eqn, dictionary()),
    dx_(dx),
    f_(0.5)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::stepUnivariateRootSolver::~stepUnivariateRootSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::stepUnivariateRootSolver::findRoot
(
    const scalar xm,
    const scalar x0,
    const scalar x1,
    const label li
) const
{
    initialise(xm);
    scalar x = x0;
    scalar dx = dx_;
    scalar yLower = eqn_.fx(x0, li);

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        scalar y = eqn_.fx(x + dx, li);
        if (y*yLower < 0)
        {
            dx *= f_;
        }
        else
        {
            x += dx;
        }
        if (converged(dx, y))
        {
            break;
        }
        printStepInformation(x);

    }

    return printFinalInformation(x);
}

// ************************************************************************* //
