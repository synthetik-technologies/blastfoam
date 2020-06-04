/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "JohnsonJacksonFrictionalStress.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace frictionalStressModels
{
    defineTypeNameAndDebug(JohnsonJackson, 0);

    addToRunTimeSelectionTable
    (
        frictionalStressModel,
        JohnsonJackson,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::
JohnsonJackson
(
    const dictionary& dict
)
:
    frictionalStressModel(dict),
    coeffDict_(dict.optionalSubDict(typeName + "Coeffs")),
    Fr_("Fr", dimensionSet(1, -1, -2, 0, 0), coeffDict_),
    eta_("eta", dimless, coeffDict_),
    p_("p", dimless, coeffDict_),
    phi_("phi", dimless, coeffDict_),
    alphaDeltaMin_("alphaDeltaMin", dimless, coeffDict_),
    alphaMinFriction_
    (
        "alphaMinFriction",
        dimless,
        coeffDict_
    )
{
    phi_ *= constant::mathematical::pi/180.0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::
~JohnsonJackson()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::
frictionalPressure
(
    const volScalarField& alphap,
    const volScalarField& alphaMax
) const
{
    return
        Fr_*pow(max(alphap - alphaMinFriction_, scalar(0)), eta_)
       /pow(max(alphaMax - alphap, alphaDeltaMin_), p_);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::
frictionalPressurePrime
(
    const volScalarField& alphap,
    const volScalarField& alphaMax
) const
{
    return Fr_*
    (
        eta_*pow(max(alphap - alphaMinFriction_, scalar(0)), eta_ - 1.0)
       *(alphaMax - alphap)
      + p_*pow(max(alphap - alphaMinFriction_, scalar(0)), eta_)
    )/pow(max(alphaMax - alphap, alphaDeltaMin_), p_ + 1.0);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::nu
(
    const phaseModel& phase,
    const volScalarField& alphap,
    const volScalarField& alphaMax,
    const volScalarField& pf,
    const volSymmTensorField& D
) const
{
    return dimensionedScalar(dimTime, 0.5)*pf*sin(phi_);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::
alphaMinFriction
(
    const volScalarField& alphap,
    const volScalarField& alphaMax
) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("alphaMinFriction", alphap.group()),
                alphap.mesh().time().timeName(),
                alphap.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            alphap.mesh(),
            alphaMinFriction_
        )
    );
}


bool Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::read()
{
    coeffDict_ <<= dict_.optionalSubDict(typeName + "Coeffs");

    Fr_.read(coeffDict_);
    eta_.read(coeffDict_);
    p_.read(coeffDict_);

    phi_.read(coeffDict_);
    phi_ *= constant::mathematical::pi/180.0;

    alphaDeltaMin_.read(coeffDict_);
    alphaMinFriction_.read(coeffDict_);

    return true;
}


// ************************************************************************* //
