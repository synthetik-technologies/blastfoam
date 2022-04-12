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

#include "FibonacciUnivariateMinimizationScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(FibonacciUnivariateMinimizationScheme, 0);
    addToRunTimeSelectionTable
    (
        minimizationScheme,
        FibonacciUnivariateMinimizationScheme,
        dictionaryUnivariate
    );
    addToRunTimeSelectionTable
    (
        univariateMinimizationScheme,
        FibonacciUnivariateMinimizationScheme,
        dictionaryZero
    );
    addToRunTimeSelectionTable
    (
        univariateMinimizationScheme,
        FibonacciUnivariateMinimizationScheme,
        dictionaryOne
    );
    addToRunTimeSelectionTable
    (
        univariateMinimizationScheme,
        FibonacciUnivariateMinimizationScheme,
        dictionaryTwo
    );
}

const Foam::scalar Foam::FibonacciUnivariateMinimizationScheme::goldenRatio =
    (sqrt(5.0) + 1.0)/2.0;

const Foam::scalar Foam::FibonacciUnivariateMinimizationScheme::s =
    (1.0 - sqrt(5.0))/(1.0 + sqrt(5.0));


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FibonacciUnivariateMinimizationScheme::FibonacciUnivariateMinimizationScheme
(
    const scalarUnivariateEquation& eqn,
    const dictionary& dict
)
:
    univariateMinimizationScheme(eqn, dict)
{
    checkY_ = true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::FibonacciUnivariateMinimizationScheme::minimize
(
    const scalar x,
    const scalar x1,
    const scalar x2,
    const label li
) const
{
    scalar a = min(x1, x2);
    scalar b = max(x1, x2);
    label n = maxSteps_;

    scalar rho = 1.0/(goldenRatio*(1.0 - pow(s, n + 1))/(1.0 - pow(s, n)));

    scalar c, yc;
    scalar d = rho*b + (1.0 - rho)*a;
    scalar yd = eqn_.fx(d, li);

    for (stepi_ = 0; stepi_ < n; stepi_++)
    {
        if (stepi_ == n - 1)
        {
            c = 0.01*a + 0.99*b;
        }
        else
        {
            c = rho*a + (1.0 - rho)*b;
        }

        eqn_.limit(c);
        yc = eqn_.fx(c, li);

        if (convergedX(a, b) && convergedY(yc, yd))
        {
            break;
        }

        if (yc < yd)
        {
            b = d;
            d = c;
            yd = yc;
        }
        else
        {
            a = b;
            b = c;
        }

        rho =
            1.0
           /(
                goldenRatio
               *(1.0 - pow(s, n + 1 - stepi_))
               /(1.0 - pow(s, n - stepi_))
            );
        printStepInformation(0.5*(a + b));
    }

    return printFinalInformation(0.5*(a + b));
}

// ************************************************************************* //
