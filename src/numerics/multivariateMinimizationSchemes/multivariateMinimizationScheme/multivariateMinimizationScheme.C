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
    const scalarList& xErrors,
    const scalarList& yErrors
) const
{
    xErrors_ = mag(xErrors);
    yErrors_ = mag(yErrors);

    bool convergedX = true;
    bool convergedY = true;
    forAll(xErrors_, i)
    {
        if (xErrors_[i] > xTolerances_[i])
        {
            convergedX = false;
            break;
        }
    }
    forAll(yErrors_, i)
    {
        if (yErrors_[i] > yTolerances_[i])
        {
            convergedY = false;
            break;
        }
    }
    return convergedX && convergedY;
}


void Foam::multivariateMinimizationScheme::printStepInformation
(
    const scalarList& vals
) const
{
    if (debug)
    {
        Info<< "Step " << stepi_ << ":" << nl
            << "    xErrors=" << xErrors_ << nl
            << "    yErrors=" << yErrors_ << nl
            << "    values: " << vals << endl;
    }
}


void Foam::multivariateMinimizationScheme::printFinalInformation() const
{
    if (stepi_ < maxSteps_ && debug)
    {
        Info<< "Converged in " << stepi_ << " iterations" << nl
            << "    Final xErrors=" << xErrors_ << nl
            << "    Final yErrors=" << yErrors_ << endl;
    }
    else if (stepi_ >= maxSteps_)
    {
        WarningInFunction
            << "Did not converge in " << maxSteps_ << " iterations" << nl 
            << "    Final xErrors= " << xErrors_ << nl
            << "    Final yErrors= " << yErrors_ << endl;
    }
}


Foam::scalar Foam::multivariateMinimizationScheme::lineSearch
(
    const scalarList& x0,
    const scalarList& grad,
    const label li,
    scalarList& fx
) const
{
    if (debug > 3)
    {
        Info<< "Conducting line search" << endl;
    }
    scalar alpha = 2.0;
    const scalarList fx0(fx);
    scalarField x1(x0 - grad);

    label iter = 0;
    do
    {
        alpha *= tau_;
        x1 = x0 - alpha*grad;
        eqns_.limit(x1);
        eqns_.f(x1, li, fx);
    } while (fx[0] > fx0[0] && iter++ < maxSteps_);
    return alpha;
}


Foam::scalar Foam::multivariateMinimizationScheme::norm(const scalarList& lst) const
{
    scalar sum = 0;
    forAll(lst, i)
    {
        sum += magSqr(lst[i]);
    }
    return sum;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multivariateMinimizationScheme::multivariateMinimizationScheme
(
    const scalarMultivariateEquation& eqns,
    const dictionary& dict
)
:
    eqns_(eqns),
    xTolerances_
    (
        dict.lookupOrDefault<scalarField>
        (
            "xTolerances",
            scalarField(eqns_.lower().size(), 1e-6)
        )
    ),
    yTolerances_
    (
        dict.lookupOrDefault<scalarField>
        (
            "yTolerances",
            scalarField(eqns_.nEqns(), 1e-6)
        )
    ),
    maxSteps_(dict.lookupOrDefault<scalar>("maxSteps", 100)),
    stepi_(0),
    xErrors_(eqns.lower().size(), great),
    yErrors_(eqns.nEqns(), great),
    tau_(dict.lookupOrDefault<scalar>("tau", 0.5))
{
    if (eqns_.nEqns() > 1)
    {
        FatalErrorInFunction
            << "Only a single output is allowed" << abort(FatalError);
    }

    if (dict.found("xTolerance"))
    {
        xTolerances_ = dict.lookup<scalar>("xTolerance");
    }
    if (dict.found("yTolerance"))
    {
        yTolerances_ = dict.lookup<scalar>("yTolerance");
    }
}


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
