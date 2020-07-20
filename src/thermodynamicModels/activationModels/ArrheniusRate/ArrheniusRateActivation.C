/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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
    activationModel(mesh, dict, phaseName),
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
      :diameterModel::New(mesh, dict, phaseName)
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
{
    if (dict.found("points") || dict.found("useCOM"))
    {
        List<vector> points(1, Zero);
        if (dict.found("useCOM"))
        {
            points[0] =
                centerOfMass
                (
                    mesh,
                    mesh.lookupObject<volScalarField>
                    (
                        IOobject::groupName("alpha", phaseName)
                    )
                );
        }
        else
        {
            points = dict.lookupType<List<vector>>("points");
        };

        scalar r(dict.lookupOrDefault("radius", 0.0));

        Info<< "Initiation Points: " << nl
            << "    " << points << endl;
        if (r > small)
        {
            Info<< "Setting cells within " << r << " m from "
                << "initiation points as activated." << endl;
        }

        forAll(points, pti)
        {
            label nCells = 0;
            if (r > small)
            {
                forAll(mesh.C(), celli)
                {
                    if (mag(mesh.C()[celli] - points[pti]) < r)
                    {
                        lambda_[celli] = 1.0;
                        nCells++;
                    }
                }
            }
            else
            {
                label celli = mesh.findCell(points[pti]);
                if (celli >= 0)
                {
                    lambda_[celli] = 1.0;
                    nCells++;
                }
            }
            if (returnReduce(nCells, sumOp<label>()) == 0)
            {
                WarningInFunction
                    << "No cells were activated using the "
                    << "detonation point " << points[pti] << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::activationModels::ArrheniusRateActivation::~ArrheniusRateActivation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::activationModels::ArrheniusRateActivation::delta() const
{
    const volScalarField& T = lambda_.mesh().lookupObject<volScalarField>(TName_);
    volScalarField R
    (
        IOobject
        (
            "R",
            lambda_.time().timeName(),
            lambda_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        lambda_.mesh(),
        dimensionedScalar("0", inv(dimTime), 0.0)
    );
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
    }
    return R*(1.0 - lambda_);
}

void Foam::activationModels::ArrheniusRateActivation::setODEFields
(
    const label nSteps,
    const labelList& oldIs,
    const label& nOld,
    const labelList& deltaIs,
    const label nDelta
)
{
    activationModel::setODEFields(nSteps, oldIs, nOld, deltaIs, nDelta);
    if (dModel_.valid())
    {
        dModel_->setODEFields(nSteps, oldIs, nOld, deltaIs, nDelta);
    }
}


void Foam::activationModels::ArrheniusRateActivation::clearODEFields()
{
    activationModel::clearODEFields();
    if (dModel_.valid())
    {
        dModel_->clearODEFields();
    }
}


void Foam::activationModels::ArrheniusRateActivation::solve
(
    const label stepi,
    const scalarList& ai,
    const scalarList& bi
)
{
    if (dModel_.valid())
    {
        dModel_->solve(stepi, ai, bi);
    }
    activationModel::solve(stepi, ai, bi);
}
// ************************************************************************* //
