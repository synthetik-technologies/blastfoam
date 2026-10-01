/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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

#include "sphericalHeatTransfer.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(sphericalHeatTransfer, 0);
    addToRunTimeSelectionTable
    (
        heatTransferModel,
        sphericalHeatTransfer,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::sphericalHeatTransfer::sphericalHeatTransfer
(
    const dictionary& dict,
    const phasePair& pair
)
:
    dispersedHeatTransferModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatTransferModels::sphericalHeatTransfer::~sphericalHeatTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::heatTransferModels::sphericalHeatTransfer::Ki() const
{
    const volScalarField& d = pair_.dispersed().d();
    const volScalarField& d0 = pair_.dispersed().d0();
    const volScalarField dEff = max(d, d0);
    return
        60.0
       *pair_.continuous().thermo().kappa()
       *d
       /pow3(dEff);
}


Foam::scalar Foam::heatTransferModels::sphericalHeatTransfer::cellKi
(
    const label celli
) const
{
    const scalar d = pair_.dispersed().celld(celli);
    const scalar d0 = pair_.dispersed().celld0(celli);
    const scalar dEff = max(d, d0);
    return
        60.0
       *pair_.continuous().thermo().kappa()[celli]
       *d
       /pow3(dEff);
}


// ************************************************************************* //
