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

#include "bisectionUnivariateMinimizationScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bisectionUnivariateMinimizationScheme, 0);
    addToRunTimeSelectionTable
    (
        univariateMinimizationScheme,
        bisectionUnivariateMinimizationScheme,
        dictionaryZero
    );
    addToRunTimeSelectionTable
    (
        univariateMinimizationScheme,
        bisectionUnivariateMinimizationScheme,
        dictionaryOne
    );
    addToRunTimeSelectionTable
    (
        univariateMinimizationScheme,
        bisectionUnivariateMinimizationScheme,
        dictionaryTwo
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bisectionUnivariateMinimizationScheme::bisectionUnivariateMinimizationScheme
(
    const equation& eqn,
    const dictionary& dict
)
:
    univariateMinimizationScheme(eqn, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::bisectionUnivariateMinimizationScheme::minimize
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
    scalar yMean= yLow;

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        if (converged(xHigh - xLow))
        {
            break;
        }

        if (yHigh < yLow)
        {
            xLow = xMean;
            yLow = eqn_.fx(xLow, li);
            yMean = yHigh;
        }
        else
        {
            xHigh = xMean;
            yHigh = eqn_.fx(xHigh, li);
            yMean = yLow;
        }

        if (converged(yMean))
        {
            break;
        }

        xMean = (xLow + xHigh)*0.5;

        printStepInformation(xMean);
    }
    converged(yMean);

    return printFinalInformation(xMean);
}

// ************************************************************************* //
