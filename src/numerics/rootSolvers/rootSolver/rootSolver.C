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

bool Foam::rootSolver::converged(const scalarField& errors) const
{
    errors_ = mag(errors);
    forAll(errors_, i)
    {
        if (errors_[i] > tolerances_[i])
        {
            return false;
        }
    }
    return true;
}


void Foam::rootSolver::printStepInformation
(
    const scalarField& vals
) const
{
    if (debug)
    {
        Info<< "Step " << stepi_
            << ", errors=" << errors_
            << ", values: " << vals << endl;
    }
}


void Foam::rootSolver::printFinalInformation() const
{
    if (stepi_ < maxSteps_ && debug)
    {
        Info<< "Converged in " << stepi_ << " iterations"
            << ", final errors=" << errors_;
    }
    else if (stepi_ >= maxSteps_)
    {
        WarningInFunction
            << "Did not converge, final errors= " << errors_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rootSolver::rootSolver
(
    const multivariateEquation<scalar>& eqns,
    const dictionary& dict
)
:
    eqns_(eqns),
    tolerances_
    (
        dict.lookupOrDefault<scalarField>
        (
            "tolerances",
            scalarField(eqns_.nVar(), 1e-6)
        )
    ),
    maxSteps_(dict.lookupOrDefault<scalar>("maxSteps", 100)),
    stepi_(0),
    errors_(eqns.nEqns(), great)
{
    if (dict.found("tolerance"))
    {
        tolerances_ = dict.lookup<scalar>("tolerance");
    }
}


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
    const scalarField& x0
) const
{
    return this->findRoots(x0, eqns_.lowerLimits(), eqns_.upperLimits(), 0);
}


Foam::tmp<Foam::scalarField> Foam::rootSolver::solve
(
    const scalarField& x0,
    const label li
) const
{
    return this->findRoots(x0, eqns_.lowerLimits(), eqns_.upperLimits(), li);
}


Foam::tmp<Foam::scalarField> Foam::rootSolver::solve
(
    const scalarField& x0,
    const scalarField& xLow,
    const scalarField& xHigh
) const
{
    return this->findRoots(x0, xLow, xHigh, 0);
}


Foam::tmp<Foam::scalarField> Foam::rootSolver::solve
(
    const scalarField& x0,
    const scalarField& xLow,
    const scalarField& xHigh,
    const label li
) const
{
    return this->findRoots(x0, xLow, xHigh, li);
}


// ************************************************************************* //
