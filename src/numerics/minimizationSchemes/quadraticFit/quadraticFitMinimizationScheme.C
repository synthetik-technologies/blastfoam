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

#include "quadraticFitMinimizationScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(quadraticFitMinimizationScheme, 0);
    addToRunTimeSelectionTable
    (
        minimizationScheme,
        quadraticFitMinimizationScheme,
        dictionaryZero
    );
    addToRunTimeSelectionTable
    (
        minimizationScheme,
        quadraticFitMinimizationScheme,
        dictionaryOne
    );
    addToRunTimeSelectionTable
    (
        minimizationScheme,
        quadraticFitMinimizationScheme,
        dictionaryTwo
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::quadraticFitMinimizationScheme::quadraticFitMinimizationScheme
(
    const scalarEquation& eqn,
    const dictionary& dict
)
:
    minimizationScheme(eqn, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::quadraticFitMinimizationScheme::solve
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

    scalar ya = eqn_.f(a, li);
    scalar yb = eqn_.f(b, li);
    scalar yc = eqn_.f(c, li);

    scalar x = b;
    scalar yx = yb;

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        x =
            0.5
           *(
                ya*(sqr(b) - sqr(c))
              + yb*(sqr(c) - sqr(a))
              + yc*(sqr(a) - sqr(b))
            )/stabilise
            (
                ya*(b - c) + yb*(c - a) + yc*(a - b),
                tolerance_
            );
        if (converged(x - b))
        {
            return x;
        }

        yx = eqn_.f(x, li);

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
    }
    return x;
}

// ************************************************************************* //
