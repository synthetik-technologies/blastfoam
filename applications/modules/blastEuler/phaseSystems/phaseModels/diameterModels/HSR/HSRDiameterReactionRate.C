/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2025
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

#include "HSRDiameterReactionRate.H"
#include "phaseModel.H"
#include "physicoChemicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterReactionRates
{
    defineTypeNameAndDebug(HSR, 0);
    addToRunTimeSelectionTable(diameterReactionRate, HSR, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterReactionRates::HSR::HSR
(
    const diameterModel& dModel,
    const dictionary& dict
)
:
    diameterReactionRate(dModel),
    A_("A", dimVelocity, dict),
    Ta_
    (
        dimTemperature,
        dict.found("Ta") || !dict.found("Ea")
      ? dict.lookup<scalar>("Ta", dimTemperature)
      : dict.lookup<scalar>("Ea", dimEnergy/dimMoles)
       /constant::physicoChemical::RR.value()
    ),
    oxidantName_(dict.lookupOrDefault("oxidant", word::null))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterReactionRates::HSR::~HSR()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::diameterReactionRates::HSR::dDdt
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    tmp<volScalarField> tmpk
    (
        volScalarField::New
        (
            typeName + ":dDdt",
            p.mesh(),
            A_
        )
    );

    const volScalarField& dp = dModel_.d();
    volScalarField Ap(dModel_.A());

    const phaseModel& phase =
        p.mesh().lookupObject<phaseModel>
        (
            IOobject::groupName
            (
                "alpha",
                dp.group()
            )
        );

    tmp<volScalarField> vDot
    (
        dModel_.A()
       *A_
       *exp(-Ta_/(phase.Ts() + dimensionedScalar(dimTemperature, small)))
    );

    if (!oxidantName_.empty())
    {
        vDot.ref() *= p.mesh().lookupObject<volScalarField>
        (
            oxidantName_
        );
    }

    return vDot/dModel_.dVdD();
}


// ************************************************************************* //
