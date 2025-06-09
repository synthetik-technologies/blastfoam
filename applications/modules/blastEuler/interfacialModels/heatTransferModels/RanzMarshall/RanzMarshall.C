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

#include "RanzMarshall.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(RanzMarshall, 0);
    addToRunTimeSelectionTable(heatTransferModel, RanzMarshall, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::RanzMarshall::RanzMarshall
(
    const dictionary& dict,
    const phasePair& pair
)
:
    NuHeatTransfer(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatTransferModels::RanzMarshall::~RanzMarshall()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::heatTransferModels::RanzMarshall::K() const
{
    return
        6.0
       *max(pair_.dispersed(), residualAlpha_)
       *pair_.continuous().thermo().kappa()
       *this->NuModel_->Nu()
       /sqr(pair_.dispersed().d());
}


Foam::scalar
Foam::heatTransferModels::RanzMarshall::cellK(const label celli) const
{
    return
        6.0
       *max(pair_.dispersed()[celli], residualAlpha_.value())
       *pair_.continuous().thermo().kappa()[celli]
       *this->NuModel_->cellNu(celli)
       /sqr(pair_.dispersed().celld(celli));
}


// ************************************************************************* //
