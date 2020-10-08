/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020 Synthetik Applied Technologies
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

#include "programmedIgnitionActivation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace activationModels
{
    defineTypeNameAndDebug(programmedIgnitionActivation, 0);
    addToRunTimeSelectionTable(activationModel, programmedIgnitionActivation, dictionary);

}
    template<>
    const char*
    NamedEnum<activationModels::programmedIgnitionActivation::burnModel, 3>::names[] =
    {
        "beta",
        "programmed",
        "programmedBeta"
    };

    const NamedEnum<activationModels::programmedIgnitionActivation::burnModel, 3>
        activationModels::programmedIgnitionActivation::burnModelNames_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::activationModels::programmedIgnitionActivation::programmedIgnitionActivation
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    activationModel(mesh, dict, phaseName),
    rho_
    (
        mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("rho", phaseName)
        )
    ),
    vDet_("vDet", dimVelocity, dict),
    Pcj_("Pcj", dimPressure, dict),
    rho0_
    (
        "rho0",
        dimDensity,
        dict.parent().subDict("products").subDict("equationOfState")
    ),
    Vcj_("Vcj", 1.0/rho0_ - Pcj_/sqr(rho0_*vDet_)),
    model_(burnModelNames_.read(dict.lookup("burnModel"))),
    detPointID_
    (
        IOobject
        (
            IOobject::groupName("detPointID", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        -1.0
    ),
    tIgn_
    (
        IOobject
        (
            IOobject::groupName("tIgn", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimTime, great)
    ),
    finished_(false)
{
    scalarField distance(detPointID_.size(), great);
    forAll(this->detonationPoints_, pointi)
    {
        const detonationPoint& dp = this->detonationPoints_[pointi];
        forAll(distance, celli)
        {
            scalar d = mag(mesh.cellCentres()[celli] - dp);
            if (d < distance[celli])
            {
                distance[celli] = d;
                detPointID_[celli] = scalar(pointi);
                tIgn_[celli] =
                    dp.delay() + mag(mesh_.C()[celli] - dp)/vDet_.value();
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::activationModels::programmedIgnitionActivation::~programmedIgnitionActivation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::activationModels::programmedIgnitionActivation::solve()
{
    if (finished_)
    {
        ddtLambda_ = volScalarField::New
        (
            "ddt(" + lambda_.name() + ")",
            lambda_.mesh(),
            dimensionedScalar(inv(dimTime), 0.0)
        );
        return;
    }

    if (min(lambda_).value() == 1.0 && this->step() == 1)
    {
        finished_ = true;
    }

    dimensionedScalar dt(this->deltaT());
    dimensionedScalar t(this->time());
    volScalarField lambdaOld(lambda_);

    // Do not include volume changes
    this->storeAndBlendOld(lambdaOld, lambdaOld_, false);

    const cellList& cells = mesh_.cells();
    const scalarField Sf(mag(mesh_.faceAreas()));
    lambda_ = lambdaOld;

    forAll(this->detonationPoints_, pointi)
    {
        detonationPoint& dp = this->detonationPoints_[pointi];
        dimensionedScalar delay(dimTime, dp.delay());
        dimensionedScalar detonationFrontDistance
        (
            max(t - delay, dimensionedScalar(dimTime, 0.0))*vDet_
        );
        dimensionedVector xDet
        (
            "xDet",
            dimLength,
            dp
        );

        forAll(lambda_, celli)
        {
            //- Activate if closest point is setActivated
            //  rounding is used for AMR
            if
            (
                mag(detPointID_[celli] - scalar(pointi)) < 1.0
             || detPointID_[celli] < 0
            )
            {
                scalar lambdaBeta = 0;
                scalar lambdaProgram = 0;

                if (model_ == BETA || model_ == PROGRAMMEDBETA)
                {
                    lambdaBeta =
                        (1.0 - 1.0/max(rho_[celli], small))
                       /(1.0 - Vcj_.value());
                }
                if (model_ == PROGRAMMED || model_ == PROGRAMMEDBETA)
                {
                    const cell& c = cells[celli];
                    scalar A = 0.0;
                    forAll(c, facei)
                    {
                        A += Sf[c[facei]];
                    }
                    scalar edgeLength = mesh_.V()[celli]/A;
                    lambdaProgram =
                        max(t.value() - tIgn_[celli], 0.0)
                       *vDet_.value()/(1.5*edgeLength);
                }

                lambda_[celli] =
                    max
                    (
                        min(max(lambdaBeta, lambdaProgram), 1.0),
                        lambda_[celli]
                    );
            }
        }
    }
    volScalarField deltaLambda(max(lambda_ - lambdaOld, 0.0)/dt);
    this->storeAndBlendDelta(deltaLambda, deltaLambda_);

    lambda_ = lambdaOld + deltaLambda*lambda_.mesh().time().deltaT();
    lambda_.maxMin(0.0, 1.0);
    lambda_.correctBoundaryConditions();

    ddtLambda_ = tmp<volScalarField>
    (
        new volScalarField
        (
            "ddt(" + lambda_.name() + ")",
            max(lambda_ - lambdaOld, 0.0)/dt
        )
    );
}

// ************************************************************************* //
