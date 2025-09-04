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
    tolerances_ = absTolerances_;
    forAll(x, i)
    {
        tolerances_[i] = max(relTolerances_[i]*mag(x[i]), tolerances_[i]);
    }
}


bool Foam::rootSolver::converged
(
    const scalarList& dx
) const
{
    bool good = true;

    forAll(dx, i)
    {
        errors_[i] = mag(dx[i]);
        if (errors_[i] > tolerances_[i])
        {
            good = false;
        }

    }
    return good;
}


bool Foam::rootSolver::converged
(
    const scalarList& x0,
    const scalarList& x1
) const
{
    bool good = true;

    forAll(x0, i)
    {
        errors_[i] = mag(x0[i] - x1[i]);
        if (errors_[i] > tolerances_[i])
        {
            good = false;
        }

    }
    return good;
}


void Foam::rootSolver::printStepInformation(const scalarList& vals) const
{
    if (debug > 2)
    {
        Info<< "Step " << stepi_
            << ", est= " << vals
            << ", error=" << errors_ << endl;
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
     && max(errors_ - tolerances_) <= 0.0;

    if (converged && debug > 1)
    {
        Info<< indent << "Converged in " << stepi_ << " iterations" << nl
            << indent << "Final errors=" << errors_ << nl
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
                << ", errors=" << errors_ <<endl;
        }
        else
        {
            WarningInFunction
                << "Did not converge in " << stepi_ << " iterations"
                << ", roots=" << vals
                << ", errors=" << errors_ << endl;
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
    boundsFac_(dict.lookupOrDefault("boundSafteyFactor", 0.5)),
    absTolerances_
    (
        dict.lookupOrDefaultBackwardsCompatible<scalarList>
        (
            {"xAbsTolerances", "absTolerances"},
            scalarList
            (
                eqns.nVar(),
                dict.lookupOrDefaultBackwardsCompatible
                (
                    {"xAbsTolerance", "absTolerance"},
                    small
                )
            )
        )
    ),
    relTolerances_
    (
        dict.lookupOrDefaultBackwardsCompatible<scalarList>
        (
            {"xTolerances", "tolerances", "relTolerances"},
            scalarList
            (
                eqns.nVar(),
                dict.lookupOrDefaultBackwardsCompatible
                (
                    {"xTolerance", "tolerance", "relTolerance"},
                    1e-6
                )
            )
        )
    ),
    tolerances_(relTolerances_),
    maxSteps_
    (
        dict.lookupOrDefaultBackwardsCompatible<label>
        (
            {"maxSteps", "maxIter"},
            100
        )
    ),
    stepi_(0),
    errors_(eqns.nVar(), great)
{}


Foam::rootSolver::rootSolver
(
    const scalarMultivariateEquation& eqns,
    const rootSolver& solver
)
:
    eqns_(eqns),
    absTolerances_(solver.absTolerances_),
    relTolerances_(solver.relTolerances_),
    tolerances_(solver.tolerances_),
    maxSteps_(solver.maxSteps_),
    stepi_(0),
    errors_(eqns.nVar(), great)
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
