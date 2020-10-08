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

#include "linearActivation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace activationModels
{
    defineTypeNameAndDebug(linearActivation, 0);
    addToRunTimeSelectionTable(activationModel, linearActivation, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::activationModels::linearActivation::linearActivation
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    activationModel(mesh, dict, phaseName),
    vDet_("vDet", dimVelocity, dict),
    finished_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::activationModels::linearActivation::~linearActivation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::activationModels::linearActivation::solve()
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
        lambda_ =
            max
            (
                pos
                (
                    detonationFrontDistance
                  - mag(lambda_.mesh().C() - xDet)
                ),
                lambda_
            );
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
