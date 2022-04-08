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

#include "HalleyUnivariateRootSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(HalleyUnivariateRootSolver, 0);
    addToRunTimeSelectionTable(univariateRootSolver, HalleyUnivariateRootSolver, dictionaryTwo);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HalleyUnivariateRootSolver::HalleyUnivariateRootSolver
(
    const scalarMultivariateEquation& eqn,
    const dictionary& dict
)
:
    univariateRootSolver(eqn, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::HalleyUnivariateRootSolver::~HalleyUnivariateRootSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::HalleyUnivariateRootSolver::findRoot
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
    scalar f = eqn_.fx(xNew, li);
    scalar fp = eqn_.dfdx(xNew, li);
    scalar fpp = eqn_.d2fdx2(xNew, li);
    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        xNew = xOld - 2.0*f*fp/stabilise(2.0*sqr(fp) - f*fpp, yTol());
        eqn_.limit(xNew);
        f = eqn_.fx(xNew, li);

        if (converged(xNew, xOld, f))
        {
            break;
        }
        fp = eqn_.dfdx(xNew, li);
        fpp = eqn_.d2fdx2(xNew, li);

        xOld = xNew;

        printStepInformation(xNew);
    }

    return printFinalInformation(xNew);
}

// ************************************************************************* //
