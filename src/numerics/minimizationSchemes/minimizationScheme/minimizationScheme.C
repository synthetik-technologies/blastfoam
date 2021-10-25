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

#include "minimizationScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(minimizationScheme, 0);
    defineRunTimeSelectionTable(minimizationScheme, dictionaryUnivariate);
    defineRunTimeSelectionTable(minimizationScheme, dictionaryMultivariate);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::minimizationScheme::converged
(
    const scalarList& errors
) const
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


void Foam::minimizationScheme::printStepInformation
(
    const scalarList& vals
) const
{
    if (debug > 1)
    {
        Info<< "Step " << stepi_ << ":" << nl
            << "    errors=" << errors_ << nl
            << "    values: " << vals << endl;
    }
}


void Foam::minimizationScheme::printFinalInformation() const
{
    if (stepi_ < maxSteps_ && debug)
    {
        Info<< "Converged in " << stepi_ << " iterations" << nl
            << "    Final errors=" << errors_ << endl;
    }
    else if (stepi_ >= maxSteps_)
    {
        WarningInFunction
            << "Did not converge in " << maxSteps_ << " iterations" << nl
            << "    Final errors= " << errors_ << endl;
    }
}


Foam::scalar Foam::minimizationScheme::lineSearch
(
    const scalarField& x0,
    const scalarField& grad,
    const label li,
    scalar& fx
) const
{
    if (debug > 3)
    {
        Info<< "Conducting line search" << endl;
    }
    scalar alpha = 2.0;
    const scalar fx0(fx);
    scalarField x1(x0 - grad);

    label iter = 0;
    do
    {
        alpha *= tau_;
        x1 = x0 - alpha*grad;
        eqns_.limit(x1);
        eqns_.f(x1, li, fx);
    } while (fx > fx0 && iter++ < maxSteps_);
    return alpha;
}


Foam::scalar Foam::minimizationScheme::norm(const scalarList& lst) const
{
    scalar sum = 0;
    forAll(lst, i)
    {
        sum += magSqr(lst[i]);
    }
    return sum;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::minimizationScheme::minimizationScheme
(
    const scalarEquation& eqns,
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
    errors_(eqns.nVar(), great),
    tau_(dict.lookupOrDefault<scalar>("tau", 0.5))
{
    if (eqns_.nEqns() > 1)
    {
        FatalErrorInFunction
            << "Only a single output is allowed" << abort(FatalError);
    }

    if (dict.found("tolerance"))
    {
        tolerances_ = dict.lookup<scalar>("tolerance");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::minimizationScheme::solve() const
{
    return solve
    (
        (eqns_.lowerLimits() + eqns_.upperLimits())*0.5,
        eqns_.lowerLimits(),
        eqns_.upperLimits(),
        0
    );
}


Foam::tmp<Foam::scalarField> Foam::minimizationScheme::solve
(
    const scalarField& x0
) const
{
    return solve(x0, eqns_.lowerLimits(), eqns_.upperLimits(), 0);
}


Foam::tmp<Foam::scalarField> Foam::minimizationScheme::solve
(
    const scalarField& x0,
    const label li
) const
{
    return this->solve(x0, eqns_.lowerLimits(), eqns_.upperLimits(), li);
}


Foam::tmp<Foam::scalarField> Foam::minimizationScheme::solve
(
    const scalarField& x0,
    const scalarField& xLow,
    const scalarField& xHigh
) const
{
    return solve(x0, xLow, xHigh, 0);
}


Foam::tmp<Foam::scalarField> Foam::minimizationScheme::solve
(
    const scalarField& x0,
    const scalarField& xLow,
    const scalarField& xHigh,
    const label li
) const
{
    return minimize(x0, xLow, xHigh, li);
}


// ************************************************************************* //
