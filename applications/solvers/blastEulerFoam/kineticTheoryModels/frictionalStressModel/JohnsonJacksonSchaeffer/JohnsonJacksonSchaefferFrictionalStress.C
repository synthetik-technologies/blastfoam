/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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

#include "JohnsonJacksonSchaefferFrictionalStress.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace frictionalStressModels
{
    defineTypeNameAndDebug(JohnsonJacksonSchaeffer, 0);

    addToRunTimeSelectionTable
    (
        frictionalStressModel,
        JohnsonJacksonSchaeffer,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::
JohnsonJacksonSchaeffer::JohnsonJacksonSchaeffer
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
    ),
    alphaMinFrictionByAlphap_
    (
        "alphaMinFrictionByAlphap",
        dimless,
        coeffDict_
    )
{
    phi_ *= constant::mathematical::pi/180.0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::
JohnsonJacksonSchaeffer::~JohnsonJacksonSchaeffer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::
JohnsonJacksonSchaeffer::frictionalPressure
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
Foam::kineticTheoryModels::frictionalStressModels::
JohnsonJacksonSchaeffer::frictionalPressurePrime
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
Foam::kineticTheoryModels::frictionalStressModels::
JohnsonJacksonSchaeffer::nu
(
    const phaseModel& phase,
    const volScalarField& alphap,
    const volScalarField& alphaMax,
    const volScalarField& pf,
    const volSymmTensorField& D
) const
{
    volScalarField alphaMinFriction(alphaMinFrictionByAlphap_*alphaMax);

    tmp<volScalarField> tnu
    (
        volScalarField::New
        (
            "JohnsonJacksonSchaeffer:nu",
            phase.mesh(),
            dimensionedScalar(dimensionSet(0, 2, -1, 0, 0), 0)
        )
    );

    volScalarField& nuf = tnu.ref();

    forAll(D, celli)
    {
        if (alphap[celli] > alphaMinFriction[celli])
        {
            nuf[celli] =
                0.5*pf[celli]*sin(phi_.value())
               /(
                    sqrt((1.0/3.0)*sqr(tr(D[celli])) - invariantII(D[celli]))
                  + SMALL
                );
        }
    }

    const fvPatchList& patches = phase.mesh().boundary();
    const volVectorField& U = phase.U();

    volScalarField::Boundary& nufBf = nuf.boundaryFieldRef();

    forAll(patches, patchi)
    {
        if (!patches[patchi].coupled())
        {
            nufBf[patchi] =
                (
                    pf.boundaryField()[patchi]*sin(phi_.value())
                   /(
                        mag(U.boundaryField()[patchi].snGrad())
                      + SMALL
                    )
                );
        }
    }

    // Correct coupled BCs
    nuf.correctBoundaryConditions();

    return tnu;
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::JohnsonJacksonSchaeffer::
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


bool Foam::kineticTheoryModels::frictionalStressModels::
JohnsonJacksonSchaeffer::read()
{
    coeffDict_ <<= dict_.optionalSubDict(typeName + "Coeffs");

    Fr_.read(coeffDict_);
    eta_.read(coeffDict_);
    p_.read(coeffDict_);

    phi_.read(coeffDict_);
    phi_ *= constant::mathematical::pi/180.0;

    alphaDeltaMin_.read(coeffDict_);

    return true;
}


// ************************************************************************* //
