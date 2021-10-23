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

#include "gradientDescentMultivariateMinimizationScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gradientDescentMultivariateMinimizationScheme, 0);
    addToRunTimeSelectionTable
    (
        multivariateMinimizationScheme,
        gradientDescentMultivariateMinimizationScheme,
        dictionaryZero
    );
    addToRunTimeSelectionTable
    (
        multivariateMinimizationScheme,
        gradientDescentMultivariateMinimizationScheme,
        dictionaryOne
    );
    addToRunTimeSelectionTable
    (
        multivariateMinimizationScheme,
        gradientDescentMultivariateMinimizationScheme,
        dictionaryTwo
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gradientDescentMultivariateMinimizationScheme::gradientDescentMultivariateMinimizationScheme
(
    const scalarMultivariateEquation& eqns,
    const dictionary& dict
)
:
    multivariateMinimizationScheme(eqns, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::gradientDescentMultivariateMinimizationScheme::minimize
(
    const scalarList& x0,
    const scalarList& xLow,
    const scalarList& xHigh,
    const label li
) const
{
    tmp<scalarField> txNew(new scalarField(x0));
    scalarField& xNew = txNew.ref();
    scalarField xOld(xNew);
    scalarField fx(eqns_.nEqns());
    eqns_.f(x0, li, fx);
    const scalarList& dx(eqns_.dx());

    scalarField grad(x0.size());
    forAll(grad, cmpti)
    {
        scalarField x1(x0);
        x1[cmpti] += dx[cmpti];
        grad[cmpti] = ((eqns_.f(x1, li) - fx)/(dx[cmpti]))()[0];
    }

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        xOld = xNew;

        scalar alpha = lineSearch(xOld, grad, li, fx);

        xNew = xOld - alpha*grad;
        eqns_.limit(xNew);

        if (converged(xNew - xOld, fx))
        {
            break;
        }

        forAll(grad, cmpti)
        {
            scalarField x1(xNew);
            x1[cmpti] += dx[cmpti];        
            grad[cmpti] = ((eqns_.f(x1, li) - fx)/(dx[cmpti]))()[0];
        }

        printStepInformation(xNew);
    }
    printFinalInformation();
    return txNew;
}

// ************************************************************************* //
