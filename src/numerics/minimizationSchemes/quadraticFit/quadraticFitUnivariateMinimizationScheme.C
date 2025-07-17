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

#include "quadraticFitUnivariateMinimizationScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace minimizationSchemes
{
namespace univariate
{
    defineTypeNameAndDebug(quadraticFit, 0);
    addToRunTimeSelectionTable
    (
        minimizationScheme,
        quadraticFit,
        dictionaryUnivariate
    );
    addToRunTimeSelectionTable
    (
        univariateMinimizationScheme,
        quadraticFit,
        dictionaryZero
    );
    addToRunTimeSelectionTable
    (
        univariateMinimizationScheme,
        quadraticFit,
        dictionaryOne
    );
    addToRunTimeSelectionTable
    (
        univariateMinimizationScheme,
        quadraticFit,
        dictionaryTwo
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::minimizationSchemes::univariate::quadraticFit::quadraticFit
(
    const scalarUnivariateEquation& eqn,
    const dictionary& dict
)
:
    univariateMinimizationScheme(eqn, dict)
{
    checkY_ = true;
}


Foam::minimizationSchemes::univariate::quadraticFit::quadraticFit
(
    const scalarUnivariateEquation& eqn,
    const quadraticFit& solver
)
:
    univariateMinimizationScheme(eqn, solver)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::minimizationSchemes::univariate::quadraticFit::minimize
(
    const scalar x0,
    const scalar x1,
    const scalar x2,
    const label li
) const
{
    scalar a = x1;
    scalar b = x0;
    scalar c = x2;

    scalar ya = eqn_.fx(a, li);
    scalar yb = eqn_.fx(b, li);
    scalar yc = eqn_.fx(c, li);

    scalar x = b;
    scalar yx = yb;

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        if
        (
            (converged(a, b, ya, yb))
         || (converged(b, c, yb, yc))
        )
        {
            break;
        }

        x =
            0.5
           *(
                ya*(sqr(b) - sqr(c))
              + yb*(sqr(c) - sqr(a))
              + yc*(sqr(a) - sqr(b))
            )/stabilise
            (
                ya*(b - c) + yb*(c - a) + yc*(a - b),
                small
            );

        eqn_.limit(x);
        yx = eqn_.fx(x, li);

        if (x > b)
        {
            if (yx > yb)
            {
                c = x;
                yc = yx;
            }
            else
            {
                a = b;
                b = x;
                ya = yb;
                yb = yx;
            }
        }
        else
        {
            if (yx > yb)
            {
                a = x;
                ya = yx;
            }
            else
            {
                c = b;
                b = x;
                yc = yb;
                yb = yx;
            }
        }

        printStepInformation(x);
    }

    return printFinalInformation(x);
}

// ************************************************************************* //
