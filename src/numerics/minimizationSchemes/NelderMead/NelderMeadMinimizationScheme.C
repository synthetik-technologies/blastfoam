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

#include "NelderMeadMinimizationScheme.H"
#include "SortableList.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(NelderMeadMinimizationScheme, 0);
    addToRunTimeSelectionTable
    (
        minimizationScheme,
        NelderMeadMinimizationScheme,
        dictionaryMultivariate
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::NelderMeadMinimizationScheme::NelderMeadMinimizationScheme
(
    const scalarEquation& eqns,
    const dictionary& dict
)
:
    minimizationScheme(eqns, dict),
    reflectionCoeff_(dict.lookupOrDefault<scalar>("reflectionCoeff", 1.0)),
    expansionCoeff_(dict.lookupOrDefault<scalar>("expansionCoeff", 2.0)),
    contractionCoeff_(dict.lookupOrDefault<scalar>("contractionCoeff", 0.5))

{
    if (reflectionCoeff_ <= 0)
    {
        FatalErrorInFunction
            << "reflectionCoeff should be greater than 0, but a value of "
            << reflectionCoeff_ << " was given" << nl
            << abort(FatalError);
    }
    if (expansionCoeff_ <= max(1.0, reflectionCoeff_))
    {
        FatalErrorInFunction
            << "expansionCoeff should be greater than "
            << max(1.0, reflectionCoeff_) << ", but a value of "
            << expansionCoeff_ << " was given" << nl
            << abort(FatalError);
    }
    if (contractionCoeff_ <= 0 || contractionCoeff_ >= 1)
    {
        FatalErrorInFunction
            << "contractionCoeff should be greater than 0 and less than 1 "
            << ", but a value of " << contractionCoeff_
            << " was given" << nl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::NelderMeadMinimizationScheme::minimize
(
    const scalarField& x0,
    const scalarField& xMin,
    const scalarField& xMax,
    const label li
) const
{
    const label np = eqns_.nVar() + 1;
    List<scalarField> points(np, x0);
    SortableList<scalar> ys(np);
    eqns_.f(points[0], li, ys[0]);
    for (label i = 1; i < np; i++)
    {
        points[i][i - 1] += (xMin[i-1] + xMax[i-1])*0.5;
        eqns_.limit(points[i]);
        eqns_.f(points[i], li, ys[i]);
    }

    // Create an indirect list to the points using the sort map
    tmp<scalarField> txMean(new scalarField(points[0]));
    scalarField& xMean = txMean.ref();
    scalarField xVar(x0.size(), 0.0);

    mean(points, xMean);
    variance(points, xMean, xVar);
    scalarField xStd(sqrt(xVar));
    if (normalize_)
    {
        xStd /= stabilise(xMean, small);
    }

    // Create some variables to reuse
    scalarField xReflection, xTmp;
    scalar yReflection, yTmp;

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        if (converged(xStd))
        {
            break;
        }

        // Sort values, point order is automatically updated since
        // the indirect list uses a reference to the sort map
        ys.sort();
        points = UIndirectList<scalarField>(points, ys.indices())();

        // Get references to the best point and its function value
        const scalarField& xLow = points[0];
        const scalar& yLow = ys[0];

        // Get references to the worst point and its function value
        scalarField& xHigh = points[np - 1];
        scalar& yHigh = ys[np - 1];

        //- Calculate the mean, neglecting the worst point
        mean(points, xMean, np-1);

        // Reflected point
        xReflection = xMean + reflectionCoeff_*(xMean - xHigh);
        if
        (
            max(pos(xReflection - tolerances_ - xMax))
         || max(neg(xReflection + tolerances_ - xMin))
        )
        {
            yReflection = great;
        }
        else
        {
            eqns_.limit(xReflection);
            eqns_.f(xReflection, li, yReflection);
        }

//         if
//         (
//             max(pos(xReflection - xMax)) > 0
//          || max(neg(xReflection - xMin)) > 0
//         )
//         {
//             Info<<xReflection<<xMax<<(xReflection - xMax)<<endl;
//             Info<<xReflection<<xMin<<(xReflection - xMin)<<endl;
//             scalarField shift
//             (
//                 max(0.0, xReflection - xMax)
//               + min(0.0, xReflection - xMin)
//             );
//             Info<<shift<<endl;
//             forAll(points, i)
//             {
//                 points[i] -= shift;
//                 eqns_.f(points[i], li, ys[i]);
//             }
//         }
        if (yReflection < yLow)
        {
            // Expansion point
            xTmp = xMean + expansionCoeff_*(xReflection - xMean);
            eqns_.limit(xTmp);
            eqns_.f(xTmp, li, yTmp);

            if (yTmp < yReflection)
            {
                xHigh = xTmp;
                yHigh = yTmp;
            }
            else
            {
                xHigh = xReflection;
                yHigh = yReflection;
            }
        }
        else if (yReflection >= ys[np-2])
        {
            if (yReflection <= yHigh)
            {
                xHigh = xReflection;
                yHigh = yReflection;
            }

            // Contraction point
            xTmp = xMean + contractionCoeff_*(xHigh - xMean);
            if
            (
                max(pos(xTmp - tolerances_ - xMax))
             || max(neg(xTmp + tolerances_ - xMin))
            )
            {
                yTmp = great;
            }
            else
            {
                eqns_.limit(xTmp);
                eqns_.f(xTmp, li, yTmp);
            }

            // Half the distance from all points to the lowest point
            if (yTmp > yHigh)
            {
                for (label i = 1; i < np; i++)
                {
                    points[i] = (points[i] + xLow)*0.5;
                    eqns_.f(points[i], li, ys[i]);
                }
            }
            else
            {
                // Set the last point and its value to the contracted point
                xHigh = xTmp;
                yHigh = yTmp;
            }
        }
        else
        {
            // Flip the points
            xHigh = xReflection;
            yHigh = yReflection;
        }

        // Update mean, variance and standard deviation
        mean(points, xMean);
        variance(points, xMean, xVar);
        xStd = sqrt(xVar);
        if (normalize_)
        {
            xStd /= stabilise(xMean, small);
        }

        printStepInformation(xMean);
    }
    xMean = points[0];
    printFinalInformation();
    return txMean;
}

// ************************************************************************* //
