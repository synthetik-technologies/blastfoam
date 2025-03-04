/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "packingLimitModel.H"
#include "SortableList.H"
#include "constantDiameterModel.H"
#include "zeroGradientFvPatchFields.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(packingLimitModel, 0);

    defineRunTimeSelectionTable(packingLimitModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::packingLimitModel::packingLimitModel
(
    const dictionary& dict,
    const masterSystem& system
)
:
    system_(system),
    mesh_(system.fluid().mesh())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::packingLimitModel::~packingLimitModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::packingLimitModel::updateAlphaMax
(
    volScalarField::Internal& alphaMaxI,
    const scalar defaultAlphaMax
) const
{
    const UPtrList<phaseModel>& phases(system_.phases());
    if (phases.size() == 1)
    {
        return;
    }

    const volScalarField& alphap = system_.alpha();
    const scalar& rAlpha = system_.residualAlpha().value();
    bool constantDiameters = true;
    forAll(phases, phasei)
    {
        if (!isA<diameterModels::constantDiameter>(phases[phasei].dModel()))
        {
            constantDiameters = false;
        }
    }

    // Only sort diameters in one cell to save time
    if (constantDiameters)
    {
        // Sort diameters from largest to smallest
        SortableList<scalar> ds(phases.size());
        forAll(phases, phasei)
        {
            ds[phasei] = phases[phasei].d()()[0];
        }
        ds.sort();

        forAll(alphaMaxI, celli)
        {
            if (alphap[celli] > rAlpha)
            {
                alphaMaxI[celli] = alphaMax(celli, ds);
            }
            else
            {
                alphaMaxI[celli] = defaultAlphaMax;
            }
        }
    }
    // Sort particle diameters for every cell
    else
    {
        SortableList<scalar> ds(phases.size());
        forAll(alphaMaxI, celli)
        {
            if (alphap[celli] > rAlpha)
            {
                forAll(phases, phasei)
                {
                    ds[phasei] = phases[phasei].celld(celli);
                }
                ds.sort();

                alphaMaxI[celli] = alphaMax(celli, ds);
            }
            else
            {
                alphaMaxI[celli] = defaultAlphaMax;
            }
        }
    }
}


// ************************************************************************* //
