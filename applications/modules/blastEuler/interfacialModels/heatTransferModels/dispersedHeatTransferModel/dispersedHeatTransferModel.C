/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-24 Jeff Heylmun:    Added return functions for acceleration
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

#include "dispersedHeatTransferModel.H"
#include "phasePair.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(dispersedHeatTransferModel, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::heatTransferModels::dispersedHeatTransferModel::dispersedHeatTransferModel
(
    const dictionary& dict,
    const phasePair& pair
)
:
    heatTransferModel(dict, pair),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        dict.lookupOrDefault<scalar>
        (
            "residualAlpha",
            pair_.ordered()
          ? pair_.dispersed().residualAlpha().value()
          : pair_.phase1().residualAlpha().value()
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatTransferModels::dispersedHeatTransferModel::~dispersedHeatTransferModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::heatTransferModels::dispersedHeatTransferModel::stabK() const
{
    return
        max
        (
            pair_.dispersed(),
            pair_.dispersed().residualAlpha()
        )*Ki();
}


Foam::scalar Foam::heatTransferModels::dispersedHeatTransferModel::cellStabK
(
    const label celli
) const
{
    return
        max
        (
            pair_.dispersed()[celli],
            pair_.dispersed().residualAlpha().value()
        )*cellKi(celli);
}


Foam::tmp<Foam::volScalarField>
Foam::heatTransferModels::dispersedHeatTransferModel::K() const
{
    return pair_.dispersed()*Ki();
}


Foam::scalar Foam::heatTransferModels::dispersedHeatTransferModel::cellK
(
    const label celli
) const
{
    return pair_.dispersed()[celli]*cellKi(celli);
}


// ************************************************************************* //
