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

#include "activationModel.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(activationModel, 0);
    defineRunTimeSelectionTable(activationModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::activationModel::activationModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    lambda_
    (
        IOobject
        (
            IOobject::groupName("lambda", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0),
        wordList(mesh.boundaryMesh().size(), "zeroGradient")
    ),
    e0_
    (
        dimensionedScalar("E0", dimPressure, dict)
       /dimensionedScalar
        (
            "rho0",
            dimDensity,
            dict.parent().subDict("products").subDict("equationOfState")
        )
    ),
    lambdaExp_(dict.lookupOrDefault("lambdaExp", 1.0)),
    alphaRhoName_
    (
        dict.lookupOrDefault
        (
            "rhoName",
            phaseName == word::null
          ? "rho"
          : IOobject::groupName("alphaRho", phaseName)
        )
    ),
    alphaRhoPhiName_
    (
        dict.lookupOrDefault
        (
            "rhoPhiName",
            phaseName == word::null
          ? "rhoPhi"
          : IOobject::groupName("alphaRhoPhi", phaseName)
        )
    ),

    maxDLambda_(dict.lookupOrDefault("maxDLambda", 1.0)),
    limit_(maxDLambda_ != 1.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::activationModel::~activationModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::activationModel::limit()
{
    if (!limit_)
    {
        return;
    }
    tmp<volScalarField> dLambda(lambda_ - lambda_.oldTime());
    lambda_ = lambda_.oldTime() + min(dLambda, maxDLambda_);
}

void Foam::activationModel::setODEFields
(
    const label nSteps,
    const labelList& oldIs,
    const label& nOld,
    const labelList& deltaIs,
    const label nDelta
)
{
    oldIs_ = oldIs;
    nOld_ = nOld;
    deltaIs_ = deltaIs;
    nDelta_ = nDelta;

    lambdaOld_.resize(nOld_);
    deltaLambda_.resize(nDelta_);
    deltaAlphaRhoLambda_.resize(nDelta_);
}


void Foam::activationModel::clearODEFields()
{
    lambdaOld_.clear();
    lambdaOld_.resize(nOld_);

    ddtLambda_.clear();

    deltaLambda_.clear();
    deltaLambda_.resize(nDelta_);

    deltaAlphaRhoLambda_.clear();
    deltaAlphaRhoLambda_.resize(nDelta_);
}


void Foam::activationModel::solve
(
    const label stepi,
    const scalarList& ai,
    const scalarList& bi
)
{
    const volScalarField& alphaRho
    (
        lambda_.mesh().lookupObject<volScalarField>(alphaRhoName_)
    );
    const surfaceScalarField& alphaRhoPhi
    (
        lambda_.mesh().lookupObject<surfaceScalarField>(alphaRhoPhiName_)
    );

    dimensionedScalar dT(alphaRho.time().deltaT());
    if (oldIs_[stepi - 1] != -1)
    {
        lambdaOld_.set
        (
            oldIs_[stepi - 1],
            new volScalarField(lambda_)
        );
    }

    volScalarField lambdaOld(lambda_*ai[stepi - 1]);
    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = oldIs_[i];
        if (fi != -1 && ai[fi] != 0)
        {
            lambdaOld += ai[fi]*lambdaOld_[fi];
        }
    }

    scalar f = bi[stepi - 1];
    for (label i = 0; i < stepi - 1; i++)
    {
        f += bi[i];
    }

    volScalarField deltaLambda(delta());
    deltaLambda = Foam::min(deltaLambda, (1.0 - lambda_)/(f*dT));
    if (deltaIs_[stepi - 1] != -1)
    {
        deltaLambda_.set
        (
            deltaIs_[stepi - 1], new volScalarField(deltaLambda)
        );
    }
    deltaLambda *= bi[stepi - 1];
    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = deltaIs_[i];
        if (fi != -1 && bi[fi] != 0)
        {
            deltaLambda += bi[fi]*deltaLambda_[fi];
        }
    }

    lambda_ = lambdaOld + deltaLambda*dT;
    lambda_.min(1);
    lambda_.max(0);
    lambda_.correctBoundaryConditions();

    if (!ddtLambda_.valid())
    {
        ddtLambda_ = tmp<volScalarField>
        (
            new volScalarField(Foam::max(lambda_ - lambdaOld, 0.0)/(f*dT))
        );
    }
    else
    {
        ddtLambda_ = Foam::max(lambda_ - lambdaOld, 0.0)/(f*dT);
    }

    volScalarField deltaAlphaRhoLambda(fvc::div(alphaRhoPhi, lambda_));
    if (deltaIs_[stepi - 1] != -1)
    {
        deltaAlphaRhoLambda_.set
        (
            deltaIs_[stepi - 1],
            new volScalarField(deltaAlphaRhoLambda)
        );
    }

    deltaAlphaRhoLambda *= bi[stepi - 1];
    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = deltaIs_[i];
        if (fi != -1 && bi[fi] != 0)
        {
            deltaAlphaRhoLambda.ref() +=
                bi[fi]*deltaAlphaRhoLambda_[fi];
        }
    }

    lambda_ =
        (
            lambdaOld*alphaRho.oldTime()
            + dT*(ddtLambda_()*f*alphaRho - deltaAlphaRhoLambda)
        )/max(alphaRho, dimensionedScalar(dimDensity, 1e-10));
    lambda_.min(1);
    lambda_.max(0);
    lambda_.correctBoundaryConditions();
}

Foam::tmp<Foam::volScalarField> Foam::activationModel::ESource() const
{
    return ddtLambda()*e0_;
}


Foam::tmp<Foam::volScalarField> Foam::activationModel::ddtLambda() const
{
    return ddtLambda_();
}

// ************************************************************************* //
