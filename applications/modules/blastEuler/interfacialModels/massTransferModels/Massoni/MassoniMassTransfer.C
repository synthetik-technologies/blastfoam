/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derived work of OpenFOAM.

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

#include "MassoniMassTransfer.H"
#include "phaseSystem.H"
#include "phasePair.H"
#include "heatTransferModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace massTransferModels
{
    defineTypeNameAndDebug(Massoni, 0);
    addToRunTimeSelectionTable(massTransferModel, Massoni, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massTransferModels::Massoni::Massoni
(
    const dictionary& dict,
    const phasePair& pair
)
:
    massTransferModel(dict, pair),
    Tsat_(saturationTemperatureModel::New("saturationTemperature", dict))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::massTransferModels::Massoni::~Massoni()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::massTransferModels::Massoni::K() const
{
    const phaseModel& phase1 = pair_.phase1();
    const phaseModel& phase2 = pair_.phase2();
    const volScalarField& p = phase1.fluid().p();
    const volScalarField Tsat(Tsat_->Tsat(p));

    const blendedHeatTransferModel& ht =
        phase1.fluid().lookupBlendedInterfacialModel<blendedHeatTransferModel>
        (
            pair_
        );
    tmp<volScalarField> tQ
    (
        volScalarField::New
        (
            "Q",
            phase1.mesh(),
            dimensionedScalar(dimDensity/dimTime*sqr(dimVelocity), 0.0)
        )
    );
    volScalarField Q = tQ.ref();
    if (ht.hasModel(phase1))
    {
        Q += ht.model(phase1).K()*(Tsat - phase1.thermo().T());
    }
    if (ht.hasModel(phase2))
    {
        Q += ht.model(phase2).K()*(Tsat - phase2.thermo().T());
    }
    return
        tQ/stabilise
        (
            phase1.thermo().ha(p, Tsat) + 0.5*magSqr(phase1.U())
          - phase2.thermo().ha(p, Tsat) - 0.5*magSqr(phase2.U()),
            dimensionedScalar(sqr(dimVelocity), small)
        );
}

// ************************************************************************* //
