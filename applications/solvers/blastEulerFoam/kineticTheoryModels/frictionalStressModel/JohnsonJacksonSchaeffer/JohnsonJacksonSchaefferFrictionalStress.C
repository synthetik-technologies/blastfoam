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
    ),
    alphaMinFrictionByAlphap_
    (
        "alphaMinFrictionByAlphap",
        dimless,
        coeffDict()
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
    const phaseModel& phase,
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
    const phaseModel& phase,
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
JohnsonJacksonSchaeffer::mu
(
    const phaseModel& phase,
    const volScalarField& alphap,
    const volScalarField& alphaMax,
    const volScalarField& pf
) const
{
    volScalarField alphaMinFriction(alphaMinFrictionByAlphap_*alphaMax);

    tmp<volScalarField> tmu
    (
        volScalarField::New
        (
            word(typeName + ":mu"),
            phase.mesh(),
            dimensionedScalar("mu", dimensionSet(1, -1, -1, 0, 0), 0.0)
        )
    );
    volScalarField& muf = tmu.ref();

    volSymmTensorField D(symm(fvc::grad(phase.U())));
    forAll(D, celli)
    {
        if (alphap[celli] > alphaMinFriction[celli])
        {
            muf[celli] =
                0.5*pf[celli]*sin(phi_.value())
               /(
                    sqrt((1.0/3.0)*sqr(tr(D[celli])) - invariantII(D[celli]))
                  + SMALL
                );
        }
    }

    const fvPatchList& patches = phase.mesh().boundary();
    const volVectorField& U = phase.U();

    volScalarField::Boundary& mufBf = muf.boundaryFieldRef();

    forAll(patches, patchi)
    {
        if (!patches[patchi].coupled())
        {
            mufBf[patchi] =
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
    muf.correctBoundaryConditions();

    return tmu;
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
    Fr_.read(coeffDict());
    eta_.read(coeffDict());
    p_.read(coeffDict());

    phi_.read(coeffDict());
    phi_ *= constant::mathematical::pi/180.0;

    alphaDeltaMin_.read(coeffDict());

    return true;
}


// ************************************************************************* //
