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
    defineRunTimeSelectionTable(minimizationScheme, dictionaryZero);
    defineRunTimeSelectionTable(minimizationScheme, dictionaryOne);
    defineRunTimeSelectionTable(minimizationScheme, dictionaryTwo);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::minimizationScheme::converged(const scalar error) const
{
    error_ = mag(error);
    return error_ < tolerance_;
}


void Foam::minimizationScheme::printStepInformation(const scalar val) const
{
    if (debug > 2)
    {
        Info<< "Step " << stepi_
            << ", error=" << error_
            << ", value: " << val << nl<< endl;
    }
}


Foam::scalar
Foam::minimizationScheme::printFinalInformation(const scalar val) const
{
    if (stepi_ < maxSteps_ && debug > 1)
    {
        Info<< "Converged in " << stepi_ << " iterations"
            << ", final error=" << error_ << endl;
    }
    else if (stepi_ >= maxSteps_ && debug)
    {
        WarningInFunction
            << "Did not converge, final error= " << error_ << endl;
    }
    return val;
}

void Foam::minimizationScheme::sample
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
        scalar y = eqn_.f(((scalar(i) + 0.5)*dx), li);
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

Foam::minimizationScheme::minimizationScheme(const scalarEquation& eqn, const dictionary& dict)
:
    eqn_(eqn),
    tolerance_(dict.lookupOrDefault<scalar>("tolerance", 1e-6)),
    maxSteps_(dict.lookupOrDefault<scalar>("maxSteps", 100)),
    nSamples_(dict.lookupOrDefault<label>("nSamples", 0)),
    stepi_(0),
    error_(great)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::minimizationScheme::solve() const
{
    return solve
    (
        (eqn_.lower() + eqn_.upper())*0.5,
        eqn_.lower(),
        eqn_.upper(),
        0
    );
}


Foam::scalar Foam::minimizationScheme::solve(const scalar x0) const
{
    return solve(x0, eqn_.lower(), eqn_.upper(), 0);
}


Foam::scalar Foam::minimizationScheme::solve
(
    const scalar x0,
    const label li
) const
{
    return this->solve(x0, eqn_.lower(), eqn_.upper(), li);
}


Foam::scalar Foam::minimizationScheme::solve
(
    const scalar x0,
    const scalar xLow,
    const scalar xHigh
) const
{
    return solve(x0, xLow, xHigh, 0);
}


Foam::scalar Foam::minimizationScheme::solve
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


void Foam::minimizationScheme::printNoConvergence() const
{
    if (debug)
    {
        WarningInFunction
            << "Did not converge with in " << stepi_ << " steps." << nl
            << "Final error=" << error_ <<endl;
    }
}


// ************************************************************************* //
