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
    const scalarList& x,
    const scalarList& xLow,
    const scalarList& xHigh,
    const label li
) const
{
    tmp<scalarField> txNew(new scalarField(x));
    RectangularMatrix<scalar> J(eqns_.nEqns(), x.size());
    scalarField& xNew = txNew.ref();
    scalarField xOld(xNew);

    scalarField fx(eqns_.nEqns());

    scalarList dx(x.size(), 1e-3);
    if (eqns_.nDerivatives() > 0)
    {
        eqns_.jacobian(x, li, fx, J);
    }
    else
    {
        J = eqns_.calculateJacobian(x, dx, li);
        fx = eqns_.f(x, li);
    }
    scalarField grad(x.size());
    for (label cmpti = 0; cmpti < x.size(); cmpti++)
    {
        grad[cmpti] = J[0][cmpti];
    }
    scalarField gradOld(grad);

    scalar alpha(0.5);
    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        xOld = xNew;
        gradOld = grad;

        xNew -= alpha*grad;
        eqns_.limit(xNew);

        if (converged(xNew - xOld))
        {
            break;
        }

        if (eqns_.nDerivatives() > 0)
        {
            eqns_.jacobian(xNew, li, fx, J);
        }
        else
        {
            J = eqns_.calculateJacobian(xNew, dx, li);
            fx = eqns_.f(xNew, li);
        }
        for (label cmpti = 0; cmpti < x.size(); cmpti++)
        {
            grad[cmpti] = J[0][cmpti];
        }
        alpha =
            mag(sum((xNew - xOld)*(grad - gradOld)))
           /max(sqrt(sum(sqr(grad - gradOld))), small);


        printStepInformation(xNew);
    }
    printFinalInformation();
    return txNew;
}

// ************************************************************************* //
