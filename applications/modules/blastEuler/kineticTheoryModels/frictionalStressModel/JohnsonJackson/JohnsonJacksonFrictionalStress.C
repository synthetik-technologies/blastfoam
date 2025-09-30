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
    const dictionary& dict,
    const masterSystem& master
)
:
    frictionalStressModel(dict, master),
    Fr_("Fr", dimensionSet(1, -1, -2, 0, 0), coeffDict()),
    eta_("eta", dimless, coeffDict()),
    p_("p", dimless, coeffDict()),
    phi_("phi", dimless, coeffDict()),
    alphaDeltaMin_("alphaDeltaMin", dimless, coeffDict()),
    alphaMinFriction_
    (
        "alphaMinFriction",
        dimless,
        coeffDict()
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
    const phaseModel& phase,
    const volScalarField& alphap,
    const volScalarField& alphaMinFriction,
    const volScalarField& alphaMax
) const
{
    return
        Fr_*pow(max(alphap - alphaMinFriction, scalar(0)), eta_)
       /pow(max(alphaMax - alphap, alphaDeltaMin_), p_);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::
frictionalPressurePrime
(
    const phaseModel& phase,
    const volScalarField& alphap,
    const volScalarField& alphaMinFriction,
    const volScalarField& alphaMax
) const
{
    volScalarField alphapMin(max(alphap - alphaMinFriction, scalar(0)));
    volScalarField alphapMax(max(alphaMax - alphap, alphaDeltaMin_));
    return
        Fr_
       *pow(alphapMin, eta_ - 1)
       /pow(alphapMax, p_ + 1)
       *(eta_*alphapMax + p_*alphapMin);

}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::mu
(
    const phaseModel& phase,
    const volScalarField& alphap,
    const volScalarField& alphaMinFriction,
    const volScalarField& alphaMax,
    const volScalarField& pf
) const
{
    return volScalarField::New
    (
        word(typeName + ":mu"),
        dimensionedScalar(dimTime, 0.5)*pf*sin(phi_)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::
alphaMinFriction
(
    const volScalarField& alphap,
    const volScalarField& alphaMax
) const
{
    return volScalarField::New
    (
        IOobject::groupName("alphaMinFriction", alphap.group()),
        alphap.mesh(),
        alphaMinFriction_
    );
}


bool Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::read()
{
    Fr_.read(coeffDict());
    eta_.read(coeffDict());
    p_.read(coeffDict());

    phi_.read(coeffDict());
    phi_ *= constant::mathematical::pi/180.0;

    alphaDeltaMin_.read(coeffDict());
    alphaMinFriction_.read(coeffDict());

    return true;
}


// ************************************************************************* //
