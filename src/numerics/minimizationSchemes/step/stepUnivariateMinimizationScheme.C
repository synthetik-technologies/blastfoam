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

#include "stepUnivariateMinimizationScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(stepUnivariateMinimizationScheme, 0);
    addToRunTimeSelectionTable
    (
        minimizationScheme,
        stepUnivariateMinimizationScheme,
        dictionaryUnivariate
    );
    addToRunTimeSelectionTable
    (
        univariateMinimizationScheme,
        stepUnivariateMinimizationScheme,
        dictionaryZero
    );
    addToRunTimeSelectionTable
    (
        univariateMinimizationScheme,
        stepUnivariateMinimizationScheme,
        dictionaryOne
    );
    addToRunTimeSelectionTable
    (
        univariateMinimizationScheme,
        stepUnivariateMinimizationScheme,
        dictionaryTwo
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::stepUnivariateMinimizationScheme::stepUnivariateMinimizationScheme
(
    const scalarUnivariateEquation& eqn,
    const dictionary& dict
)
:
    univariateMinimizationScheme(eqn, dict),
    dx_
    (
        dict.lookupOrDefault<scalar>
        (
            "dx",
            (eqn_.upper() - eqn_.lower())/100.0
        )
    ),
    f_(dict.lookupOrDefault<scalar>("f", 0.5))
{
    checkY_ = true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::stepUnivariateMinimizationScheme::minimize
(
    const scalar,
    const scalar x0,
    const scalar x1,
    const label li
) const
{
    scalar xLower = x0;
    scalar dx = dx_;
    scalar xMid = x0 + dx*f_;
    scalar xUpper = x0 + dx;

    scalar yLower = eqn_.fx(xLower, li);
    scalar yMid = eqn_.fx(xMid, li);
    scalar yUpper = eqn_.fx(xUpper, li);

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        if (yUpper < yLower)
        {
            if (yMid > yUpper)
            {
                xLower = xUpper;
                yLower = yUpper;
            }
            else
            {
                xLower = xMid;
                yLower = yMid;
                dx *= f_;
            }
        }
        else
        {
            dx *= f_;
        }

        xMid = xLower + dx*f_;
        yMid = eqn_.fx(xMid, li);

        xUpper = xLower + dx;

        eqn_.limit(xUpper);
        yUpper = eqn_.fx(xUpper, li);

        if (convergedX(xLower, xUpper) && convergedY(yLower, yUpper))
        {
            break;
        }
        printStepInformation(xUpper);
    }

    convergedY(yUpper, yLower);
    return printFinalInformation(xLower);
}

// ************************************************************************* //
