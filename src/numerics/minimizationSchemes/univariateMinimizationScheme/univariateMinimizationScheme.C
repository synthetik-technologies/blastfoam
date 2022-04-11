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

#include "univariateMinimizationScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(univariateMinimizationScheme);
    defineRunTimeSelectionTable(univariateMinimizationScheme, dictionaryZero);
    defineRunTimeSelectionTable(univariateMinimizationScheme, dictionaryOne);
    defineRunTimeSelectionTable(univariateMinimizationScheme, dictionaryTwo);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::univariateMinimizationScheme::convergedXScale
(
    const scalar error,
    const scalar s
) const
{
    xErrors_[0] = mag(error);
    xRelErrors_[0] = xErrors_[0]/max(mag(s), small);
    return converged(xErrors_, xRelErrors_, xTolerances_, xRelTolerances_);
}


bool Foam::univariateMinimizationScheme::convergedX
(
    const scalar x1,
    const scalar x2
) const
{
    xErrors_[0] = mag(x2 - x1);
    xRelErrors_[0] = xErrors_[0]/stabilise(min(mag(x1), mag(x2)), small);
    return converged(xErrors_, xRelErrors_, xTolerances_, xRelTolerances_);
}



bool Foam::univariateMinimizationScheme::convergedYScale
(
    const scalar error,
    const scalar s
) const
{
    yErrors_[0] = mag(error);
    yRelErrors_[0] = yErrors_[0]/max(mag(s), small);
    return converged(yErrors_, yRelErrors_, yTolerances_, yRelTolerances_);
}


bool Foam::univariateMinimizationScheme::convergedY
(
    const scalar y1,
    const scalar y2
) const
{
    yErrors_[0] = mag(y2 - y1);
    yRelErrors_[0] = yErrors_[0]/stabilise(min(mag(y1), mag(y2)), small);
    return converged(yErrors_, yRelErrors_, yTolerances_, yRelTolerances_);
}


void Foam::univariateMinimizationScheme::printStepInformation
(
    const scalar val
) const
{
    if (debug > 2)
    {
        DebugInfo<< "Step: " << stepi_ << ":" << nl
            << "    Error (abs/rel): "
            << xErrors_[0] << ", " << xRelErrors_[0] << endl;
        if (checkY_)
        {
            Info<< "    Delta (abs/rel): "
                << yErrors_[0] << ", " << yRelErrors_[0] << endl;
        }
        Info<< "    Minimum: " << val << endl;
    }
}


Foam::scalar
Foam::univariateMinimizationScheme::printFinalInformation(const scalar val) const
{
    if (!debug)
    {
        return val;
    }
    bool converged =
        (xErrors_[0] - xTolerances_[0] <= 0.0)
     || (xRelErrors_[0] - xRelTolerances_[0] <= 0.0)
     || (yErrors_[0] - yTolerances_[0] <= 0.0)
     || (yRelErrors_[0] - yRelTolerances_[0] <= 0.0);

    if (converged)
    {
        Info<< "Converged in " << stepi_ << " iterations" << endl;
    }
    else if (stepi_ < maxSteps_)
    {
        WarningInFunction
            << "Did not converge due to bounds, tried "
            << stepi_ << " iterations" << endl;
    }
    else
    {
        WarningInFunction
            << "Did not converge in "
            << stepi_ << " iterations" << endl;
    }
    Info<< "    Final error (abs/rel): "
        << xErrors_[0] << ", " << xRelErrors_[0] << endl;
    if (checkY_)
    {
        Info<< "    Final delta (abs/rel): "
            << yErrors_[0] << ", " << yRelErrors_[0] << endl;
    }
    Info<< "    Minimum: " << val << endl;
    return val;
}

void Foam::univariateMinimizationScheme::sample
(
    scalar& x0,
    scalar& x1,
    const label li
) const
{
    if (nSample_ < 1)
    {
        return;
    }
    if (debug > 1)
    {
        Info<<"Pre sampling interval" << endl;
    }

    scalar dx = (x1 - x0)/scalar(nSample_ + 1);
    label minI = 0;
    scalar minY = great;
    for (label i = 0; i < nSample_; i++)
    {
        scalar y = eqn_.fx(((scalar(i) + 0.5)*dx + x0), li);
        if (y < minY)
        {
            minY = y;
            minI = i;
        }
    }
    x0 += scalar(minI)*dx;
    x1 = x0 + dx;

    if (debug > 1)
    {
        Info<<"Found minimum values in (" << x0 << "," << x1 << ")"
            << " minY = " << minY
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::univariateMinimizationScheme::univariateMinimizationScheme
(
    const scalarUnivariateEquation& eqn,
    const dictionary& dict
)
:
    minimizationScheme(eqn, dict),
    eqn_(dynamicCast<const scalarEquation>(eqn)),
    nSample_(dict.lookupOrDefault<label>("nSample", 0))
{
    nSamples_ = 0;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::univariateMinimizationScheme::solve(const scalar x0) const
{
    return solve(x0, eqn_.lower(), eqn_.upper(), 0);
}


Foam::scalar Foam::univariateMinimizationScheme::solve
(
    const scalar x0,
    const label li
) const
{
    return this->solve(x0, eqn_.lower(), eqn_.upper(), li);
}


Foam::scalar Foam::univariateMinimizationScheme::solve
(
    const scalar x0,
    const scalar xLow,
    const scalar xHigh
) const
{
    return solve(x0, xLow, xHigh, 0);
}


Foam::scalar Foam::univariateMinimizationScheme::solve
(
    const scalar x0,
    const scalar xLow,
    const scalar xHigh,
    const label li
) const
{
    scalar xMin = xLow;
    scalar xMax = xHigh;
    sample(xMin, xMax, li);
    return minimize(x0, xMin, xMax, li);
}


Foam::tmp<Foam::scalarField> Foam::univariateMinimizationScheme::minimize
(
    const scalarList& x0,
    const scalarList& xLow,
    const scalarList& xHigh,
    const label li
) const
{
    return tmp<scalarField>
    (
        new scalarField(1, solve(x0[0], xLow[0], xHigh[0], li))
    );
}

// ************************************************************************* //
