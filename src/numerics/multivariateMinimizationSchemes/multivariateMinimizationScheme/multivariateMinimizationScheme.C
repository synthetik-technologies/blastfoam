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

#include "multivariateMinimizationScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multivariateMinimizationScheme, 0);
    defineRunTimeSelectionTable(multivariateMinimizationScheme, dictionaryZero);
    defineRunTimeSelectionTable(multivariateMinimizationScheme, dictionaryOne);
    defineRunTimeSelectionTable(multivariateMinimizationScheme, dictionaryTwo);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::multivariateMinimizationScheme::converged
(
    const scalarList& errors
) const
{
    errors_ = mag(errors);
    error_ = max(errors_);

    forAll(errors_, i)
    {
        if (errors_[i] > tolerances_[i])
        {
            return false;
        }
    }
    return true;
}


void Foam::multivariateMinimizationScheme::printStepInformation
(
    const scalarList& vals
) const
{
    if (debug)
    {
        Info<< "Step " << stepi_
            << ", errors=" << errors_
            << ", values: " << vals << endl;
    }
}


void Foam::multivariateMinimizationScheme::printFinalInformation() const
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

Foam::multivariateMinimizationScheme::multivariateMinimizationScheme
(
    const scalarMultivariateEquation& eqns,
    const dictionary& dict
)
:
    eqns_(eqns),
    tolerances_
    (
        dict.lookupOrDefault<scalarField>
        (
            "tolerances",
            scalarField(eqns_.nEqns(), 1e-6)
        )
    ),
    maxSteps_(dict.lookupOrDefault<scalar>("maxSteps", 100)),
    stepi_(0),
    errors_(eqns.nEqns(), great),
    error_(great)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::multivariateMinimizationScheme::solve() const
{
    return solve
    (
        (eqns_.lower() + eqns_.upper())*0.5,
        eqns_.lower(),
        eqns_.upper(),
        0
    );
}


Foam::tmp<Foam::scalarField> Foam::multivariateMinimizationScheme::solve
(
    const scalarList& x0
) const
{
    return solve(x0, eqns_.lower(), eqns_.upper(), 0);
}


Foam::tmp<Foam::scalarField> Foam::multivariateMinimizationScheme::solve
(
    const scalarList& x0,
    const label li
) const
{
    return this->solve(x0, eqns_.lower(), eqns_.upper(), li);
}


Foam::tmp<Foam::scalarField> Foam::multivariateMinimizationScheme::solve
(
    const scalarList& x0,
    const scalarList& xLow,
    const scalarList& xHigh
) const
{
    return solve(x0, xLow, xHigh, 0);
}


Foam::tmp<Foam::scalarField> Foam::multivariateMinimizationScheme::solve
(
    const scalarList& x0,
    const scalarList& xLow,
    const scalarList& xHigh,
    const label li
) const
{
    return minimize(x0, xLow, xHigh, li);
}


// ************************************************************************* //
