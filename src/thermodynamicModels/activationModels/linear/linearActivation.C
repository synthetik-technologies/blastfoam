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
    vDet_("vDet", dimVelocity, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::activationModels::linearActivation::~linearActivation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::activationModels::linearActivation::clearODEFields()
{
    ddtLambda_.clear();
}


void Foam::activationModels::linearActivation::solve
(
    const label stepi,
    const scalarList& ai,
    const scalarList& bi
)
{
    if (stepi > 1)
    {
        return;
    }

    dimensionedScalar dt(lambda_.time().deltaT());
    volScalarField lambdaOld(lambda_);

    forAll(this->detonationPoints_, pointi)
    {
        detonationPoint& dp = this->detonationPoints_[pointi];
        dp.setActivated(lambda_, true);
        dimensionedScalar delay(dimTime, dp.delay());
        dimensionedScalar detonationFrontDistance
        (
            (lambda_.time() + dt - delay)*vDet_
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
                pos(detonationFrontDistance - mag(lambda_.mesh().C() - xDet)),
                lambda_
            );

    }
    lambda_ = max(lambdaOld, lambda_);

    ddtLambda_ = tmp<volScalarField>
    (
        new volScalarField((lambda_ - lambdaOld)/dt)
    );
}

// ************************************************************************* //
