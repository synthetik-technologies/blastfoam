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

template<class T>
class listMaxEqOp
{
public:

    void operator()(List<T>& x, const List<T>& y) const
    {
        forAll(x, i)
        {
            x[i] = max(x[i], y[i]);
        }
    }
};
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::packingLimitModel::packingLimitModel
(
    const dictionary& dict,
    const masterSystem& system
)
:
    system_(system),
    mesh_(system.fluid().mesh()),
    minAlphaMax_(1.0)
{
    forAll(system_.phases(), phasei)
    {
        minAlphaMax_ = min(minAlphaMax_, system_.phases()[phasei].alphaMax());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::packingLimitModel::~packingLimitModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::packingLimitModel::updateAlphaMax
(
    volScalarField::Internal& alphaMaxI
) const
{
    const UPtrList<phaseModel>& phases(system_.phases());
    if (phases.size() == 1)
    {
        alphaMaxI = minAlphaMax_;
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
    SortableList<scalar> ds(phases.size(), 0.0);
    if (constantDiameters)
    {
        // Sort diameters from largest to smallest
        forAll(phases, phasei)
        {
            if (alphap.size() && Pstream::myProcNo() == 1)
            {
                ds[phasei] = phases[phasei].celld(0);
            }
        }
        ds.reverseSort();

        Pstream::combineGather(ds, listMaxEqOp<scalar>());
        Pstream::scatter(ds);
    }

    forAll(alphaMaxI, celli)
    {
        if (alphap[celli] > rAlpha)
        {
            forAll(phases, phasei)
            {
                ds[phasei] = phases[phasei].celld(celli);
            }
            ds.sort();

            const scalar aM = alphaMax(celli, ds, constantDiameters);
            alphaMaxI[celli] = aM < 0 ? minAlphaMax_ : aM;
        }
        else
        {
            alphaMaxI[celli] = minAlphaMax_;
        }
    }
}


// ************************************************************************* //
