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

#include "cyclicCoordinateMinimizationScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace minimizationSchemes
{
    defineTypeNameAndDebug(cyclicCoordinate, 0);
    addToRunTimeSelectionTable
    (
        minimizationScheme,
        cyclicCoordinate,
        dictionaryMultivariate
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::minimizationSchemes::cyclicCoordinate::cyclicCoordinate
(
    const scalarUnivariateEquation& eqns,
    const dictionary& dict
)
:
    basis(eqns, dict),
    accelerate_(dict.lookupOrDefault<bool>("accelerate", true))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::minimizationSchemes::cyclicCoordinate::minimize
(
    const scalarList& x0,
    const scalarList& xMin,
    const scalarList& xMax,
    const label li
) const
{
    tmp<scalarField> txNew(new scalarField(x0));
    scalarField& xNew = txNew.ref();
    scalarField xOld(xNew);
    scalarField delta(xNew);

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        xOld = xNew;
        forAll(xNew, diri)
        {
            basis::searchDir(xNew, li, diri, 1, xNew);
        }

        delta = xNew - xOld;
        if (accelerate_)
        {
            lineSearcher().search(xNew, delta, li, xNew);
            delta = xNew - xOld;
        }
        if (convergedXScale(delta, xNew))
        {
            break;
        }
        printStepInformation(xNew);
    }
    printFinalInformation(xNew);
    return txNew;
}

// ************************************************************************* //
