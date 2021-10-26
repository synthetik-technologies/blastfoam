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

#include "NewtonRaphsonUnivariateRootSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(NewtonRaphsonUnivariateRootSolver, 0);
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        NewtonRaphsonUnivariateRootSolver,
        dictionaryOne
    );
    addToRunTimeSelectionTable
    (
        univariateRootSolver,
        NewtonRaphsonUnivariateRootSolver,
        dictionaryTwo
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::NewtonRaphsonUnivariateRootSolver::NewtonRaphsonUnivariateRootSolver
(
    const multivariateEquation<scalar>& eqn,
    const dictionary& dict
)
:
    univariateRootSolver(eqn, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::NewtonRaphsonUnivariateRootSolver::~NewtonRaphsonUnivariateRootSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::NewtonRaphsonUnivariateRootSolver::findRoot
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
        xNew = xOld - eqn_.fx(xOld, li)/stabilise(eqn_.dfdx(xOld, li), small);
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
