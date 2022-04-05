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

#include "SteffensenUnivariateRootSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SteffensenUnivariateRootSolver, 0);
    addToRunTimeSelectionTable
    (
        rootSolver,
        SteffensenUnivariateRootSolver,
        dictionaryUnivariate
    );
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        SteffensenUnivariateRootSolver,
        dictionaryZero
    );
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        SteffensenUnivariateRootSolver,
        dictionaryOne
    );
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        SteffensenUnivariateRootSolver,
        dictionaryTwo
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SteffensenUnivariateRootSolver::SteffensenUnivariateRootSolver
(
    const scalarMultivariateEquation& eqn,
    const dictionary& dict
)
:
    univariateRootSolver(eqn, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SteffensenUnivariateRootSolver::~SteffensenUnivariateRootSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::SteffensenUnivariateRootSolver::findRoot
(
    const scalar x0,
    const scalar x1,
    const scalar x2,
    const label li
) const
{
    scalar xOld = x0;
    scalar xNew = x0;
    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        scalar fx = eqn_.fx(xOld, li);
        scalar gx = eqn_.fx(xOld + fx, li)/stabilise(fx, small) - 1.0;
        xNew = xOld - eqn_.fx(xOld, li)/stabilise(gx, small);
        eqn_.limit(xNew);

        if (converged(xNew - xOld))
        {
            break;
        }
        xOld = xNew;

        printStepInformation(xNew);

    }

    return printFinalInformation(xNew);
}

// ************************************************************************* //
