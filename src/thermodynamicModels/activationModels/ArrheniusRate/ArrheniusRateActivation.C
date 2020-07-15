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
    dp_(diameterModel::New(mesh, dict, phaseName)),
    Ts_("Ts", dimTemperature, dict),
    ALow_("ALow", inv(sqr(dimLength)*dimTime), dict),
    TaLow_("TaLow", dimTemperature, dict),
    AHigh_("AHigh", inv(dimTime), dict),
    TaHigh_("TaHigh", dimTemperature, dict),
    T_
    (
        mesh.foundObject<volScalarField>(IOobject::groupName("T", phaseName))
      ? mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("T", phaseName)
        )
      : mesh.lookupObject<volScalarField>("T")
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

    forAll(R, celli)
    {
        if (T_[celli] > Ts_.value())
        {
            R[celli] =
                sqr(dp_->d()[celli])
               *AHigh_.value()
               *exp(-TaHigh_.value()/T_[celli]);
        }
        else
        {
            R[celli] =
                ALow_.value()
               *exp(-TaLow_.value()/T_[celli]);
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
    dp_->setODEFields(nSteps, oldIs, nOld, deltaIs, nDelta);
}


void Foam::activationModels::ArrheniusRateActivation::clearODEFields()
{
    activationModel::clearODEFields();
    dp_->clearODEFields();
}


void Foam::activationModels::ArrheniusRateActivation::solve
(
    const label stepi,
    const scalarList& ai,
    const scalarList& bi
)
{
    dp_->solve(stepi, ai, bi);
    activationModel::solve(stepi, ai, bi);
}
// ************************************************************************* //
