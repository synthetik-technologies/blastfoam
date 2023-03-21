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

#include "HookeJeevesMinimizationScheme.H"
#include "SortableList.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace minimizationSchemes
{
    defineTypeNameAndDebug(HookeJeeves, 0);
    addToRunTimeSelectionTable
    (
        minimizationScheme,
        HookeJeeves,
        dictionaryMultivariate
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::minimizationSchemes::HookeJeeves::HookeJeeves
(
    const scalarUnivariateEquation& eqns,
    const dictionary& dict
)
:
    minimizationScheme(eqns, dict),
    gamma_(dict.lookupOrDefault<scalar>("gamma", 0.9))

{

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::minimizationSchemes::HookeJeeves::minimize
(
    const scalarList& x0,
    const scalarList& xMin,
    const scalarList& xMax,
    const label li
) const
{
    tmp<scalarField> txNew(new scalarField(x0));
    scalarField& xNew = txNew.ref();
    scalarField xOld(x0);
    scalarField xBest(x0);
    scalar yBest = eqns_.fX(x0, li);
    scalarField alpha(xMax);
    forAll(alpha, i)
    {
        alpha[i] = (xMax[i] - xMin[i])/2.0;
    };

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        bool improved = false;
        xOld = xNew;

        forAll(xOld, diri)
        {
            xNew = xOld;
            xNew[diri] += alpha[diri];
            eqns_.limit(xNew);
            scalar y = eqns_.fX(xNew,  li);
            if (y < yBest)
            {
                xBest = xNew;
                yBest = y;
                improved = true;
            }

            xNew[diri] -= 2.0*alpha[diri];
            eqns_.limit(xNew);
            y = eqns_.fX(xNew, li);
            if (y < yBest)
            {
                xBest = xNew;
                yBest = y;
                improved = true;
            }
        }

        if (!improved)
        {
            alpha *= gamma_;
        }
        else
        {
            xNew = xBest;
            if (convergedXScale(xNew - xOld, xNew))
            {
                break;
            }
        }

        printStepInformation(xNew);
    }

    xNew.transfer(xBest);

    printFinalInformation(xNew);
    return txNew;
}

// ************************************************************************* //
