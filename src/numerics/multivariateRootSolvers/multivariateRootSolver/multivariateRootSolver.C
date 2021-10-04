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

void Foam::multivariateRootSolver::printNoConvergence() const
{
    if (debug)
    {
        WarningInFunction
            << "Did not converge with in " << stepi_ << " steps." << nl
            << "Final errors=" << errors_ <<endl;
    }
}

// ************************************************************************* //
