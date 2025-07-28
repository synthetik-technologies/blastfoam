/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
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

#include "bisectionUnivariateMinimizationScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace minimizationSchemes
{
namespace univariate
{
    defineTypeNameAndDebug(bisection, 0);
    addToRunTimeSelectionTable
    (
        minimizationScheme,
        bisection,
        dictionaryUnivariate
    );
    addToRunTimeSelectionTable
    (
        univariateMinimizationScheme,
        bisection,
        dictionaryZero
    );
    addToRunTimeSelectionTable
    (
        univariateMinimizationScheme,
        bisection,
        dictionaryOne
    );
    addToRunTimeSelectionTable
    (
        univariateMinimizationScheme,
        bisection,
        dictionaryTwo
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::minimizationSchemes::univariate::bisection::bisection
(
    const scalarUnivariateEquation& eqn,
    const dictionary& dict
)
:
    univariateMinimizationScheme(eqn, dict)
{
    checkY_ = true;
}


Foam::minimizationSchemes::univariate::bisection::bisection
(
    const scalarUnivariateEquation& eqn,
    const bisection& solver
)
:
    univariateMinimizationScheme(eqn, solver)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::minimizationSchemes::univariate::bisection::minimize
(
    const scalar x,
    const scalar x1,
    const scalar x2,
    const label li
) const
{
    scalar xLow = x1;
    scalar xHigh = x2;
    scalar xMean = 0.5*(x1 + x2);
    scalar yLow = eqn_.fx(xLow, li);
    scalar yHigh = eqn_.fx(xHigh, li);

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        if (converged(xLow, xHigh, yLow, yHigh))
        {
            break;
        }

        if (yHigh < yLow)
        {
            xLow = xMean;
            yLow = eqn_.fx(xLow, li);
        }
        else
        {
            xHigh = xMean;
            yHigh = eqn_.fx(xHigh, li);
        }

        xMean = (xLow + xHigh)*0.5;

        printStepInformation(xMean);
    }

    return printFinalInformation(xMean);
}

// ************************************************************************* //
