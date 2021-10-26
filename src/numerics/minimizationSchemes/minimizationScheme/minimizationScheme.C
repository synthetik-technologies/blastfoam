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
        Info<< "Step: " << stepi_ << ":" << nl
            << "    " << errorName() << "s: " << errors_ << nl
            << "    values: " << vals << endl;
    }
}


void Foam::minimizationScheme::printFinalInformation() const
{
    if (stepi_ < maxSteps_ && debug)
    {
        Info<< "Converged in " << stepi_ << " iterations" << nl
            << "    Final " << errorName() << "s=" << errors_ << endl;
    }
    else if (stepi_ >= maxSteps_)
    {
        WarningInFunction
            << "Did not converge in " << maxSteps_ << " iterations" << nl
            << "    Final " << errorName() << "s: " << errors_ << endl;
    }
}


Foam::scalar Foam::minimizationScheme::alpha
(
    const scalarField& grad,
    const scalarField& gradOld
) const
{
    return 1.0;
}


Foam::scalar Foam::minimizationScheme::lineSearch
(
    const scalarField& x0,
    const scalarField& grad,
    const scalarField& gradOld,
    const label li,
    scalar& fx
) const
{
    if (debug > 3)
    {
        Info<< "Conducting line search" << endl;
    }

//     scalar alpha = norm(grad)/stabilise(norm(gradOld), small);
//     scalar alpha =
//         max
//         (
//             inner(grad, grad - gradOld)/stabilise(norm(gradOld), small),
//             0.0
//         );
    scalar a = alpha(grad, gradOld);
    const scalar fx0(fx);
    scalarField x1(x0 - a*grad);
    eqns_.limit(x1);
    eqns_.f(x1, li, fx);

    label iter = 0;
    while (fx >= fx0 && iter++ < maxSteps_)
    {
        a *= tau_;
        x1 = x0 - a*grad;
        eqns_.limit(x1);
        eqns_.f(x1, li, fx);
    }
    return a;
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
    scalar alpha = 10.0;

    const scalar fx0(fx);
    scalarField x1(x0 - alpha*grad);
    eqns_.limit(x1);
    eqns_.f(x1, li, fx);

    label iter = 0;
    while (fx >= fx0 && iter++ < maxSteps_)
    {
        alpha *= tau_;
        x1 = x0 - alpha*grad;
        eqns_.limit(x1);
        eqns_.f(x1, li, fx);
    }
    return alpha;
}


Foam::scalar Foam::minimizationScheme::norm(const scalarList& lst) const
{
    scalar sum = 0;
    forAll(lst, i)
    {
        sum += magSqr(lst[i]);
    }
    return sqrt(sum/scalar(lst.size()));
}


Foam::scalar Foam::minimizationScheme::inner
(
    const scalarList& lst1,
    const scalarList& lst2
) const
{
    scalar sum = 0;
    forAll(lst1, i)
    {
        sum += lst1[i]*lst2[i];
    }
    return sum;
}


void Foam::minimizationScheme::sample
(
    scalarField& xLow,
    scalarField& xHigh,
    labelList& xBest,
    scalar& yBest,
    labelList& indicies,
    const label li,
    const label diri
) const
{
    if ((debug || minimizationScheme::debug) && diri == 0)
    {
        Info<<"Pre sampling interval" << endl;
    }

    if (diri >= xBest.size())
    {
        scalar y;
        eqns_.f((xLow + xHigh)*0.5, li, y);

        if (y < yBest)
        {
            yBest = y;
            xBest = indicies;
        }
        return;
    }


    scalar dx = (xHigh[diri] - xLow[diri])/scalar(nSamples_[diri] + 1);
    for (label i = 0; i < nSamples_[diri]; i++)
    {
        scalarField x0(xLow);
        x0[diri] += dx*scalar(i);
        scalarField x1(xHigh);
        x1[diri] += dx;
        indicies[diri] = i;
        sample(x0, x1, xBest, yBest, indicies, li, diri + 1);
    }

    if (diri != 0)
    {
        return;
    }

    forAll(xLow, i)
    {
        scalar dx = (xHigh[i] - xLow[i])/scalar(nSamples_[i] + 1);
        xLow[i] += scalar(xBest[i])*dx;
        xHigh[i] = xLow[i] + dx;
    }

    if (debug || minimizationScheme::debug)
    {
        Info<< "Found minimum values in (" << xLow << "," << xHigh << ")"
            << " yBest: " << yBest << endl;
    }
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
    nSamples_
    (
        dict.lookupOrDefault<labelList>
        (
            "nSamples",
            labelList(eqns.nVar(), 0)
        )
    ),
    normalize_(dict.lookupOrDefault("normalizeError", true)),
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
    scalarField x1(xLow);
    scalarField x2(xHigh);
    scalarField xStart(x0);
    if (min(nSamples_) > 0)
    {
        labelList xBest(x1.size(), 0);
        labelList indicies(xBest);
        scalar yBest = great;
        sample(x1, x2, xBest, yBest, indicies, li, 0);
        xStart = 0.5*(x1 + x2);
    }
    return minimize(xStart, x1, x2, li);
}

// ************************************************************************* //
