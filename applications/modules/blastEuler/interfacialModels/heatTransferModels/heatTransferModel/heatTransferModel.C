/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2025-06-09 Jeff Heylmun:    Added cell based returns
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

#include "heatTransferModel.H"
#include "BlendedInterfacialModel.H"
#include "phasePair.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heatTransferModel, 0);
    defineBlendedInterfacialModelTypeNameAndDebug(heatTransferModel, 0);
    defineRunTimeSelectionTable(heatTransferModel, dictionary);
}

const Foam::dimensionSet Foam::heatTransferModel::dimK(1, -1, -3, -1, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModel::heatTransferModel
(
    const dictionary& dict,
    const phasePair& pair
)
:
    pair_(pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatTransferModel::~heatTransferModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::blendedHeatTransferModel::stabK() const
{
    return this->evaluate
    (
        &heatTransferModel::stabK,
        "K",
        heatTransferModel::dimK,
        false
    );
}


Foam::scalar Foam::blendedHeatTransferModel::cellStabK(const label celli) const
{
    return this->evaluate(&heatTransferModel::cellStabK, false, celli);
}


Foam::tmp<Foam::volScalarField> Foam::blendedHeatTransferModel::K() const
{
    return this->evaluate
    (
        &heatTransferModel::K,
        "K",
        heatTransferModel::dimK,
        false
    );
}


Foam::scalar Foam::blendedHeatTransferModel::cellK(const label celli) const
{
    return this->evaluate(&heatTransferModel::cellK, false, celli);
}


// ************************************************************************* //
