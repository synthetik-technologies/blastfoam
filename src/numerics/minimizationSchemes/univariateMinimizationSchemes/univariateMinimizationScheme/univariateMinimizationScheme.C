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
    defineTypeNameAndDebug(univariateMinimizationScheme, 0);
    defineRunTimeSelectionTable(univariateMinimizationScheme, dictionaryZero);
    defineRunTimeSelectionTable(univariateMinimizationScheme, dictionaryOne);
    defineRunTimeSelectionTable(univariateMinimizationScheme, dictionaryTwo);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::univariateMinimizationScheme::converged(const scalar error) const
{
    errors_[0] = mag(error);
    return errors_[0] < tolerances_[0];
}


void Foam::univariateMinimizationScheme::printStepInformation(const scalar val) const
{
    if (debug > 2)
    {
        Info<< "Step " << stepi_
            << ", error=" << errors_[0]
            << ", value: " << val << nl<< endl;
    }
}


Foam::scalar
Foam::univariateMinimizationScheme::printFinalInformation(const scalar val) const
{
    if (stepi_ < maxSteps_ && debug > 1)
    {
        Info<< "Converged in " << stepi_ << " iterations"
            << ", final error=" << errors_[0] << endl;
    }
    else if (stepi_ >= maxSteps_ && debug)
    {
        WarningInFunction
            << "Did not converge, final error= " << errors_[0] << endl;
    }
    return val;
}

void Foam::univariateMinimizationScheme::sample
(
    scalar& x0,
    scalar& x1,
    const label li
) const
{
    if (nSamples_ < 1)
    {
        return;
    }
    if (debug)
    {
        Info<<"Pre sampling interval" << endl;
    }

    scalar dx = (x1 - x0)/scalar(nSamples_ + 1);
    label minI = 0;
    scalar minY = great;
    for (label i = 0; i < nSamples_; i++)
    {
        scalar y = eqn_.fx(((scalar(i) + 0.5)*dx), li);
        if (y < minY)
        {
            minY = y;
            minI = i;
        }
    }
    x0 = scalar(minI)*dx;
    x1 = x0 + dx;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::univariateMinimizationScheme::univariateMinimizationScheme
(
    const scalarEquation& eqn,
    const dictionary& dict
)
:
    minimizationScheme(eqn, dict),
    eqn_(dynamicCast<const equation>(eqn)),
    nSamples_(dict.lookupOrDefault<label>("nSamples", 0))
{}

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
    const scalarField& x0,
    const scalarField& xLow,
    const scalarField& xHigh,
    const label li
) const
{
    return tmp<scalarField>
    (
        new scalarField(1, solve(x0[0], xLow[0], xHigh[0], li))
    );
}

// ************************************************************************* //
