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

#include "multivariateRootSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multivariateRootSolver, 0);
    defineRunTimeSelectionTable(multivariateRootSolver, dictionaryZero);
    defineRunTimeSelectionTable(multivariateRootSolver, dictionaryOne);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::multivariateRootSolver::converged(const scalarField& errors) const
{
    errors_ = mag(errors);
    error_ = max(errors_);
    return error_ < tolerance_;
}


void Foam::multivariateRootSolver::printStepInformation
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


void Foam::multivariateRootSolver::printFinalInformation() const
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

Foam::multivariateRootSolver::multivariateRootSolver
(
    const scalarMultivariateEquation& eqns,
    const dictionary& dict
)
:
    eqns_(eqns),
    tolerance_(dict.lookupOrDefault<scalar>("tolerance", 1e-6)),
    maxSteps_(dict.lookupOrDefault<scalar>("maxSteps", 100)),
    stepi_(0),
    errors_(eqns.nEqns(), great),
    error_(great)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::multivariateRootSolver::solve() const
{
    return this->findRoots
    (
        ((eqns_.lower() + eqns_.upper())*0.5)(),
        eqns_.lower(),
        eqns_.upper(),
        0
    );
}


Foam::tmp<Foam::scalarField> Foam::multivariateRootSolver::solve
(
    const scalarList& x0
) const
{
    return this->findRoots(x0, eqns_.lower(), eqns_.upper(), 0);
}


Foam::tmp<Foam::scalarField> Foam::multivariateRootSolver::solve
(
    const scalarList& x0,
    const label li
) const
{
    return this->findRoots(x0, eqns_.lower(), eqns_.upper(), li);
}


Foam::tmp<Foam::scalarField> Foam::multivariateRootSolver::solve
(
    const scalarList& x0,
    const scalarList& xLow,
    const scalarList& xHigh
) const
{
    return this->findRoots(x0, xLow, xHigh, 0);
}


Foam::tmp<Foam::scalarField> Foam::multivariateRootSolver::solve
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
