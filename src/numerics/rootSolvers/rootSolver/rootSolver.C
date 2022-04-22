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

#include "rootSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rootSolver, 0);
    defineRunTimeSelectionTable(rootSolver, dictionaryUnivariate);
    defineRunTimeSelectionTable(rootSolver, dictionaryZero);
    defineRunTimeSelectionTable(rootSolver, dictionaryOne);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::rootSolver::initialise(const scalarList& x) const
{
    forAll(x, i)
    {
        xRelTols_[i] = max(xTols_[i]*mag(x[i]), xAbsTols_[i]);
    }
}


bool Foam::rootSolver::converged
(
    const scalarList& dx,
    const scalarList& y
) const
{
    bool good = true;

    forAll(y, i)
    {
        yErrors_[i] = mag(y[i]);
        if (yErrors_[i] > yTols_[i] )
        {
           good = false;
        }

    }

    forAll(dx, i)
    {
        xErrors_[i] = mag(dx[i]);
        if (xErrors_[i] > xRelTols_[i])
        {
            good = false;
        }

    }
    return good;
}


bool Foam::rootSolver::converged
(
    const scalarList& x0,
    const scalarList& x1,
    const scalarList& y
) const
{
    bool good = true;

    forAll(y, i)
    {
        yErrors_[i] = mag(y[i]);
        if (yErrors_[i] > yTols_[i] )
        {
           good = false;
        }

    }

    bool zeroDiff = true;
    forAll(x0, i)
    {
        xErrors_[i] = mag(x0[i] - x1[i]);
        if (xErrors_[i] > xRelTols_[i])
        {
            good = false;
            if (xErrors_[i] > vSmall)
            {
                zeroDiff = false;
            }
        }

    }
    return zeroDiff || good;
}


void Foam::rootSolver::printStepInformation(const scalarList& vals) const
{
    if (debug > 2)
    {
        Info<< "Step " << stepi_
            << ", est= " << vals
            << ", error=" << xErrors_ << "/" << yErrors_ << endl;
    }
}

void Foam::rootSolver::printFinalInformation(const scalarList& vals) const
{
    if (!debug)
    {
        return;
    }

    bool converged =
        (stepi_ < maxSteps_)
     && max(xErrors_ - xRelTols_) <= 0.0
     && max(yErrors_ - yTols_) <= 0.0;

    if (converged && debug > 1)
    {
        Info<< indent << "Converged in " << stepi_ << " iterations" << nl
            << indent << "Final x errors=" << xErrors_ << nl
            << indent << "Final y errors=" << yErrors_ << nl
            << indent << "Roots=" << vals << endl;
    }
    else if (!converged)
    {
        if (stepi_ < maxSteps_)
        {
            WarningInFunction
                << "Did not converge due to bounds"
                << ", tried " << stepi_ << " iterations"
                << ", est=" << vals
                << ", errors=" << xErrors_ << "/" << yErrors_ << endl;
        }
        else
        {
            WarningInFunction
                << "Did not converge in " << stepi_ << " iterations"
                << ", roots=" << vals
                << ", errors=" << xErrors_ << "/" << yErrors_ << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rootSolver::rootSolver
(
    const scalarMultivariateEquation& eqns,
    const dictionary& dict
)
:
    eqns_(eqns),
    xTols_
    (
        dict.lookupOrDefault<scalarList>
        (
            "xTolerances",
            scalarList
            (
                eqns.nVar(),
                dict.lookupOrDefault("xTolerance", 1e-6)
            )
        )
    ),
    yTols_
    (
        dict.lookupOrDefault<scalarList>
        (
            "yTolerances",
            scalarList
            (
                eqns.nEqns(),
                dict.lookupOrDefault("yTolerance", 1e-6)
            )
        )
    ),
    xAbsTols_
    (
        dict.lookupOrDefault<scalarList>
        (
            "xAbsTolerances",
            scalarList
            (
                eqns.nVar(),
                dict.lookupOrDefault("xAbsTolerance", 1e-6)
            )
        )
    ),
    xRelTols_(xTols_),
    maxSteps_(dict.lookupOrDefault<scalar>("maxSteps", 100)),
    stepi_(0),
    xErrors_(eqns.nVar(), great),
    yErrors_(eqns.nEqns(), great)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rootSolver::~rootSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::rootSolver::solve() const
{
    return this->findRoots
    (
        ((eqns_.lowerLimits() + eqns_.upperLimits())*0.5)(),
        eqns_.lowerLimits(),
        eqns_.upperLimits(),
        0
    );
}


Foam::tmp<Foam::scalarField> Foam::rootSolver::solve
(
    const scalarList& x0
) const
{
    return this->findRoots(x0, eqns_.lowerLimits(), eqns_.upperLimits(), 0);
}


Foam::tmp<Foam::scalarField> Foam::rootSolver::solve
(
    const scalarList& x0,
    const label li
) const
{
    return this->findRoots(x0, eqns_.lowerLimits(), eqns_.upperLimits(), li);
}


Foam::tmp<Foam::scalarField> Foam::rootSolver::solve
(
    const scalarList& x0,
    const scalarList& xLow,
    const scalarList& xHigh
) const
{
    return this->findRoots(x0, xLow, xHigh, 0);
}


Foam::tmp<Foam::scalarField> Foam::rootSolver::solve
(
    const scalarList& x0,
    const scalarList& xLow,
    const scalarList& xHigh,
    const label li
) const
{
    return this->findRoots(x0, xLow, xHigh, li);
}


// ************************************************************************* //
