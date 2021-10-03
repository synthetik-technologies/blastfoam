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

#include "stepRootSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(stepRootSolver, 0);
    addToRunTimeSelectionTable
    (
        rootSolver,
        stepRootSolver,
        dictionaryZero
    );
    addToRunTimeSelectionTable
    (
        rootSolver,
        stepRootSolver,
        dictionaryOne
    );
    addToRunTimeSelectionTable
    (
        rootSolver,
        stepRootSolver,
        dictionaryTwo
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::stepRootSolver::stepRootSolver
(
    const scalarEquation& eqn,
    const dictionary& dict
)
:
    rootSolver(eqn, dict),
    dx_(dict.lookupOrDefault<scalar>("dx", (eqn_.upper() - eqn_.lower())/100.0))
{}


Foam::stepRootSolver::stepRootSolver
(
    const scalarEquation& eqn,
    const scalar dx
)
:
    rootSolver(eqn, dictionary()),
    dx_(dx)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::stepRootSolver::solve
(
    const scalar x0,
    const scalar x1,
    const scalar x2,
    const label li
) const
{
    scalar x = x1;
    scalar dx = dx_;
    scalar yLower = eqn_.f(x0, li);

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        scalar y = eqn_.f(x + dx, li);
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
            return x;
        }

    }
    WarningInFunction
        << "Could not converge to the given root." << endl;

    return x;
}

// ************************************************************************* //
