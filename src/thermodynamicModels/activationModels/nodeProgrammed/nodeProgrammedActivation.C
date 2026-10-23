/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2025
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

#include "nodeProgrammedActivation.H"
#include "addToRunTimeSelectionTable.H"
#include "timeIntegrationSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace activationModels
{
    defineTypeNameAndDebug(nodeProgrammedActivation, 0);
    addToRunTimeSelectionTable(activationModel, nodeProgrammedActivation, dictionary);
}
    template<>
    const char*
    NamedEnum<activationModels::nodeProgrammedActivation::burnModel, 2>::names[] =
    {
        "programmed",
        "programmedBeta"
    };

    const NamedEnum<activationModels::nodeProgrammedActivation::burnModel, 2>
        activationModels::nodeProgrammedActivation::burnModelNames_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::activationModels::nodeProgrammedActivation::nodeProgrammedActivation
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    activationModel
    (
        mesh,
        dict,
        phaseName,
        burnModelNames_[dict.lookupOrDefault("burnModel", word("programmed"))] == PROGRAMMED
      ? 1
      : 0
    ),
    rho_
    (
        mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("rho", phaseName)
        )
    ),
    pName_(dict.lookupOrDefault("pName", phaseName != word::null ? IOobject::groupName("p", phaseName) : word("p"))),
    p_
    (
        mesh.foundObject<volScalarField>(pName_)
      ? mesh.lookupObject<volScalarField>(pName_)
      : mesh.lookupObject<volScalarField>("p")
    ),
    pMin_("pMin", dimPressure, dict.lookupOrDefault<scalar>("pMin", 1.5e5)),
    alphaPtr_
    (
        phaseName != word::null && mesh.foundObject<volScalarField>(IOobject::groupName("alpha", phaseName))
      ? &mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("alpha", phaseName)
        )
      : nullptr
    ),
    residualAlpha_(dict.lookupOrDefault<scalar>("residualAlpha", 1e-4)),
    vDet_("vDet", dimVelocity, dict),
    Pcj_("Pcj", dimPressure, dict),
    rho0_
    (
        "rho0",
        dimDensity,
        dict.parent().subDict("products").subDict("equationOfState")
    ),
    Vcj_("Vcj", 1.0 - Pcj_/(rho0_*sqr(vDet_))),
    model_(burnModelNames_[dict.lookupOrDefault("burnModel", word("programmed"))]),
    expansionBurnRatio_(dict.lookupOrDefault<scalar>("expansionBurnRatio", 0)),
    tStart_
    (
        IOobject
        (
            IOobject::groupName("tStart", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimTime, great)
    ),
    tEnd_
    (
        IOobject
        (
            IOobject::groupName("tEnd", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimTime, great)
    ),
    useDelay_(dict.lookupOrDefault("delayOffset", false)),
    timeIndex_(-1)
{
    calcArrivalTimes();
}


// * * * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * * //

void Foam::activationModels::nodeProgrammedActivation::calcArrivalTimes()
{
    const fvMesh& mesh = lambda_.mesh();
    const pointField& pts = mesh.points();
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();

    tStart_ = dimensionedScalar(dimTime, great);
    tEnd_ = dimensionedScalar(dimTime, great);

    forAll(this->detonationPoints_, pointi)
    {
        const detonationPoint& dp = this->detonationPoints_[pointi];
        const scalar delay = (useDelay_ ? dp.delay() : 0.0);

        forAll(cells, celli)
        {
            const labelList cellPts = cells[celli].labels(faces);
            scalar tMinCell = great;
            scalar tMaxCell = 0.0;

            forAll(cellPts, pi)
            {
                const label ptIndex = cellPts[pi];
                const scalar tp = delay + mag(pts[ptIndex] - dp)/vDet_.value();
                tMinCell = min(tMinCell, tp);
                tMaxCell = max(tMaxCell, tp);
            }

            if (tMinCell < tStart_[celli])
            {
                tStart_[celli] = dp.delay() + tMinCell;
                tEnd_[celli] = dp.delay() + tMaxCell;
            }
        }
    }

    timeIndex_ = mesh.time().timeIndex();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::activationModels::nodeProgrammedActivation::~nodeProgrammedActivation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::activationModels::nodeProgrammedActivation::delta() const
{
    return volScalarField::New
    (
        IOobject::groupName(type() + ":R", lambda_.group()),
        lambda_.mesh(),
        dimensionedScalar
        (
            "delta",
            inv(dimTime),
            0.0
        )
    );
}

void Foam::activationModels::nodeProgrammedActivation::correct
(
    volScalarField& lambda
)
{
    const fvMesh& mesh = lambda.mesh();
    if (mesh.topoChanged() && mesh.time().timeIndex() != timeIndex_)
    {
        calcArrivalTimes();
    }

    dimensionedScalar t(timeIntegrationSystem::t());

    forAll(lambda, celli)
    {
        const scalar alphai = alphaPtr_ ? (*alphaPtr_)[celli] : 1.0;

        if (alphai <= residualAlpha_)
        {
            lambda[celli] = 0.0;
            continue;
        }

        const scalar tStart = tStart_[celli];
        const scalar tEnd = tEnd_[celli];
        const scalar dtCell = max(tEnd - tStart, 1e-15);

        scalar lambdaProg = 0.0;
        if (t.value() >= tEnd)
        {
            lambdaProg = 1.0;
        }
        else if (t.value() > tStart)
        {
            lambdaProg = (t.value() - tStart)/dtCell;
        }

        if (model_ == PROGRAMMEDBETA)
        {
            if (lambdaProg > 0.0 && rho_[celli] > rho0_.value())
            {
                lambda[celli] = max
                (
                    max(lambda[celli], lambdaProg),
                    (1.0 - rho0_.value()/rho_[celli])/(1.0 - Vcj_.value())
                );
            }
            else
            {
                lambda[celli] = max(lambda[celli], lambdaProg);
            }
        }
        else
        {
            lambda[celli] = max(lambda[celli], lambdaProg);
        }

        if (rho_[celli] < expansionBurnRatio_*rho0_.value())
        {
            lambda[celli] = 1.0;
        }
    }
}

// ************************************************************************* //
