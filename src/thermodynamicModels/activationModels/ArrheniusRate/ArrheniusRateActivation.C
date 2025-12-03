/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "ArrheniusRateActivation.H"
#include "thermodynamicConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace activationModels
{
    defineTypeNameAndDebug(ArrheniusRateActivation, 0);
    addToRunTimeSelectionTable(activationModel, ArrheniusRateActivation, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::activationModels::ArrheniusRateActivation::ArrheniusRateActivation
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    activationModel(mesh, dict, phaseName, 0),
    rho0_
    (
        "rho0",
        dimDensity,
        dict.parent().subDict("products").subDict("equationOfState")
    ),
    dModel_
    (
        mesh.foundObject<volScalarField>(IOobject::groupName("d", phaseName))
      ? autoPtr<diameterModel>()
      : diameterModel::New(mesh, dict, phaseName)
    ),
    dp_(mesh.lookupObject<volScalarField>(IOobject::groupName("d", phaseName))),
    Tign_("Tign", dimTemperature, dict),
    Ts_("Ts", dimTemperature, dict),
    ALow_("ALow", inv(sqr(dimLength)*dimTime), dict),
    EaLow_("EaLow", dimEnergy/dimMass, dict),
    AHigh_("AHigh", inv(dimTime), dict),
    EaHigh_("EaHigh", dimEnergy/dimMass, dict),
    TName_
    (
        dict.lookupOrDefault
        (
            "TName",
            mesh.foundObject<volScalarField>(IOobject::groupName("T", phaseName))
          ? IOobject::groupName("T", phaseName)
          : "T"
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::activationModels::ArrheniusRateActivation::~ArrheniusRateActivation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::activationModels::ArrheniusRateActivation::delta() const
{
    const volScalarField& T = lambda_.mesh().lookupObject<volScalarField>(TName_);
    tmp<volScalarField> tR
    (
        volScalarField::New
        (
            IOobject::groupName(type() + ":R", lambda_.group()),
            lambda_.mesh(),
            dimensionedScalar("0", inv(dimTime), 0.0)
        )
    );
    volScalarField& R = tR.ref();
    scalar specieR(Foam::constant::thermodynamic::RR);

    forAll(R, celli)
    {
        if (T[celli] < Tign_.value())
        {
            R[celli] = 0.0;
        }
        else if (T[celli] > Ts_.value())
        {
            R[celli] =
                sqr(dp_[celli])
               *AHigh_.value()
               *exp(-EaHigh_.value()/(specieR*T[celli]));
        }
        else
        {
            R[celli] =
                ALow_.value()
               *exp(-EaLow_.value()/(specieR*T[celli]));
        }
        R[celli] *= (1.0 - lambda_[celli]);
    }

    return tR;
}


void Foam::activationModels::ArrheniusRateActivation::update()
{
    activationModel::update();
    if (dModel_.valid())
    {
        dModel_->update();
    }
}


void Foam::activationModels::ArrheniusRateActivation::solve()
{
    if (dModel_.valid())
    {
        dModel_->solve();
    }
    activationModel::solve();
}


void Foam::activationModels::ArrheniusRateActivation::solveExplicit()
{
    if (dModel_.valid())
    {
        dModel_->solveExplicit();
    }
    activationModel::solveExplicit();
}


void Foam::activationModels::ArrheniusRateActivation::storeExplicit()
{
    if (dModel_.valid())
    {
        dModel_->storeExplicit();
    }
    activationModel::storeExplicit();
}


void Foam::activationModels::ArrheniusRateActivation::solveImplicit()
{
    if (dModel_.valid())
    {
        dModel_->solveImplicit();
    }
    activationModel::solveImplicit();
}


void Foam::activationModels::ArrheniusRateActivation::clear()
{
    if (dModel_.valid())
    {
        dModel_->clear();
    }
    activationModel::clear();
}

// ************************************************************************* //
