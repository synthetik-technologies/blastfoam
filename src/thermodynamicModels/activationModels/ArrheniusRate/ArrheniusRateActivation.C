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
    rho0_("rho0", dimDensity, dict.parent()),
    dp_("dp", dimLength, dict),
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
{}


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
                sqr(dp_.value())
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

// ************************************************************************* //
