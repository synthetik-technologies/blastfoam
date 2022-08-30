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

#include "gradientDescentMinimizationScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gradientDescentMinimizationScheme, 0);
    addToRunTimeSelectionTable
    (
        minimizationScheme,
        gradientDescentMinimizationScheme,
        dictionaryUnivariate
    );
    addToRunTimeSelectionTable
    (
        minimizationScheme,
        gradientDescentMinimizationScheme,
        dictionaryMultivariate
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gradientDescentMinimizationScheme::gradientDescentMinimizationScheme
(
    const scalarUnivariateEquation& eqns,
    const dictionary& dict
)
:
    minimizationScheme(eqns, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::gradientDescentMinimizationScheme::minimize
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
    scalar fx(eqns_.fX(x0, li));
    scalarField grad(x0.size(), 0.0);
    eqns_.dfdX(x0, li, grad);
    scalarField gradOld(grad);
    scalar alpha;
    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        xOld = xNew;
        if (stepi_ <= 1)
        {
        alpha = lineSearch(xOld, grad, li, fx);
        }
        else
        {
            alpha = lineSearch(xOld, grad, gradOld, li, fx);
        }
        xNew = xOld - alpha*grad;

        eqns_.limit(xNew);

        if (convergedXScale(xNew, xOld))
        {
            break;
        }

        gradOld = grad;
        fx = eqns_.fX(xNew, li);
        eqns_.dfdX(xNew, li, grad);
        printStepInformation(xNew);
    }
    printFinalInformation(xNew);
    return txNew;
}

// ************************************************************************* //
