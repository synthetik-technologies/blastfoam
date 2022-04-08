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

#include "univariateRootSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(univariateRootSolver);
    defineRunTimeSelectionTable(univariateRootSolver, dictionaryZero);
    defineRunTimeSelectionTable(univariateRootSolver, dictionaryOne);
    defineRunTimeSelectionTable(univariateRootSolver, dictionaryTwo);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::univariateRootSolver::initialise(const scalar x) const
{
    xRelTols_[0] = max(xTols_[0]*mag(x), xAbsTols_[0]);
}


bool Foam::univariateRootSolver::converged
(
    const scalar dx,
    const scalar y
) const
{
    xErrors_[0] = mag(dx);
    yErrors_[0] = mag(y);
    return xErrors_[0] < xTols_[0] && yErrors_[0] < yTols_[0];
}


bool Foam::univariateRootSolver::converged
(
    const scalar x0,
    const scalar x1,
    const scalar y
) const
{
    xErrors_[0] = mag(x1 - x0);
    yErrors_[0] = mag(y);
    return xErrors_[0] < xRelTols_[0] && yErrors_[0] < yTols_[0];
}



void Foam::univariateRootSolver::printStepInformation(const scalar val) const
{
    if (debug > 2)
    {
        Info<< "Step " << stepi_
            << ", est= " << val
            << ", error=" << xErrors_[0] << "/" << yErrors_[0] << endl;
    }
}


Foam::scalar
Foam::univariateRootSolver::printFinalInformation(const scalar val) const
{
    if (!debug)
    {
        return val;
    }

    bool converged =
        (stepi_ < maxSteps_)
     && xErrors_[0] - xRelTols_[0] <= 0.0
     && yErrors_[0] - yTols_[0] <= 0.0;

    if (converged && debug > 1)
    {
        Info<< indent << "Converged in " << stepi_ << " iterations" << nl
            << indent << "Final x error=" << xErrors_[0] << nl
            << indent << "Final y error=" << yErrors_[0] << nl
            << indent << "Root=" << val << endl;
    }
    else if (!converged)
    {
        if (stepi_ < maxSteps_)
        {
            WarningInFunction
                << "Did not converge due to bounds"
                << ", tried " << stepi_ << " iterations"
                << ", est=" << val
                << ", error=" << xErrors_[0] << "/" << yErrors_[0] << endl;
        }
        else
        {
            WarningInFunction
                << "Did not converge in " << stepi_ << " iterations"
                << ", roots=" << val
                << ", errors=" << xErrors_[0] << "/" << yErrors_[0] << endl;
        }
    }
    return val;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::univariateRootSolver::univariateRootSolver
(
    const scalarMultivariateEquation& eqn,
    const dictionary& dict
)
:
    rootSolver(eqn, dict),
    eqn_(dynamicCast<const scalarEquation>(eqn))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::univariateRootSolver::~univariateRootSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::univariateRootSolver::solve(const scalar x0) const
{
    return solve(x0, eqn_.lower(), eqn_.upper(), 0);
}


Foam::scalar Foam::univariateRootSolver::solve
(
    const scalar x0,
    const label li
) const
{
    return solve(x0, eqn_.lower(), eqn_.upper(), li);
}


Foam::scalar Foam::univariateRootSolver::solve
(
    const scalar x0,
    const scalar xLow,
    const scalar xHigh
) const
{
    return solve(x0, xLow, xHigh, 0);
}


Foam::scalar Foam::univariateRootSolver::solve
(
    const scalar x0,
    const scalar xLow,
    const scalar xHigh,
    const label li
) const
{
    return this->findRoot(x0, xLow, xHigh, li);
}


Foam::tmp<Foam::scalarField> Foam::univariateRootSolver::findRoots
(
    const scalarField& x0,
    const scalarField& xLow,
    const scalarField& xHigh,
    const label li
) const
{
    return tmp<scalarField>
    (
        new scalarField(1, this->findRoot(x0[0], xLow[0], xHigh[0], li))
    );
}


// ************************************************************************* //
