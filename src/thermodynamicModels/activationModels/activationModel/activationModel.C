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
    integrationSystem
    (
        IOobject::groupName("activationModel", phaseName),
        mesh
    ),
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
        dict.found("E0")
      ? dimensionedScalar("E0", dimPressure, dict)
       /dimensionedScalar
        (
            "rho0",
            dimDensity,
            dict.parent().subDict("products").subDict("equationOfState")
        )
      : dimensionedScalar("e0", dimEnergy/dimMass, dict)
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
{
    this->lookupAndInitialize();
    Info<<lambdaOld_.size()<<endl;
    Info<<deltaLambda_.size()<<endl;
    Info<<deltaAlphaRhoLambda_.size()<<endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::activationModel::~activationModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::activationModel::centerOfMass
(
    const fvMesh& mesh,
    const volScalarField& alpha
) const
{
    scalarField Vtot(mesh.V()*alpha.primitiveField());
    vectorField m1(Vtot*mesh.C().primitiveField());
    return gSum(m1)/gSum(Vtot);
}


void Foam::activationModel::limit()
{
    if (!limit_)
    {
        return;
    }
    tmp<volScalarField> dLambda(lambda_ - lambda_.oldTime());
    lambda_ = lambda_.oldTime() + min(dLambda, maxDLambda_);
}


void Foam::activationModel::clearODEFields()
{
    this->clearOld(lambdaOld_);
    this->clearDelta(deltaLambda_);
    this->clearDelta(deltaAlphaRhoLambda_);
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
    volScalarField lambdaOld(lambda_);
    this->storeOld(stepi, lambdaOld, lambdaOld_);
    this->blendOld(stepi, lambdaOld, lambdaOld_, ai);


    scalar f = this->f(stepi, bi);

    volScalarField deltaLambda(delta());
    deltaLambda = Foam::min(deltaLambda, (1.0 - lambda_)/(f*dT));
    this->storeDelta(stepi, deltaLambda, deltaLambda_);
    this->blendDelta(stepi, deltaLambda, deltaLambda_, bi);

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
    this->storeDelta(stepi, deltaAlphaRhoLambda, deltaAlphaRhoLambda_);
    this->blendDelta(stepi, deltaAlphaRhoLambda, deltaAlphaRhoLambda_, bi);

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
