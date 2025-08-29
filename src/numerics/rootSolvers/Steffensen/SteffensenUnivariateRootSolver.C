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

#include "SteffensenUnivariateRootSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace rootSolvers
{
namespace univariate
{
    defineTypeNameAndDebug(Steffensen, 0);
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        Steffensen,
        dictionaryZero
    );
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        Steffensen,
        dictionaryOne
    );
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        Steffensen,
        dictionaryTwo
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rootSolvers::univariate::Steffensen::Steffensen
(
    const scalarMultivariateEquation& eqn,
    const dictionary& dict
)
:
    univariateRootSolver(eqn, dict)
{}


Foam::rootSolvers::univariate::Steffensen::Steffensen
(
    const scalarMultivariateEquation& eqn,
    const Steffensen& solver
)
:
    univariateRootSolver(eqn, solver)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rootSolvers::univariate::Steffensen::~Steffensen()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::rootSolvers::univariate::Steffensen::findRoot
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

    scalar fx = eqn_.fx(xOld, li);
    scalar gx = eqn_.fx(xOld + fx, li)/stabilise(fx, small) - 1.0;

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        xNew = xOld - eqn_.fx(xOld, li)/stabilise(gx, small);
        eqn_.limit(xNew);

        fx = eqn_.fx(xNew, li);
        if (converged(xNew, xOld))
        {
            break;
        }

        gx = eqn_.fx(xNew + fx, li)/stabilise(fx, small) - 1.0;

        xOld = xNew;
        printStepInformation(xNew);

    }

    return printFinalInformation(xNew);
}

// ************************************************************************* //
