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

#include "SchaefferFrictionalStress.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace frictionalStressModels
{
    defineTypeNameAndDebug(Schaeffer, 0);

    addToRunTimeSelectionTable
    (
        frictionalStressModel,
        Schaeffer,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::Schaeffer::Schaeffer
(
    const dictionary& dict,
    const kineticTheorySystem& kt
)
:
    frictionalStressModel(dict, kt),
    phi_("phi", dimless, coeffDict()),
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

Foam::kineticTheoryModels::frictionalStressModels::Schaeffer::~Schaeffer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::Schaeffer::
frictionalPressure
(
    const phaseModel& phase,
    const volScalarField& alphap,
    const volScalarField& alphaMax
) const
{
    volScalarField alphaMinFriction(alphaMinFrictionByAlphap_*alphaMax);

    return
        dimensionedScalar(dimensionSet(1, -1, -2, 0, 0), 1e24)
       *pow(Foam::max(alphap - alphaMinFriction, scalar(0)), 10.0);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::Schaeffer::
frictionalPressurePrime
(
    const phaseModel& phase,
    const volScalarField& alphap,
    const volScalarField& alphaMax
) const
{
    volScalarField alphaMinFriction(alphaMinFrictionByAlphap_*alphaMax);

    return
        dimensionedScalar(dimensionSet(1, -1, -2, 0, 0), 1e25)
       *pow(Foam::max(alphap - alphaMinFriction, scalar(0)), 9.0);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::Schaeffer::mu
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
Foam::kineticTheoryModels::frictionalStressModels::Schaeffer::
alphaMinFriction
(
    const volScalarField& alphap,
    const volScalarField& alphaMax
) const
{
    return alphaMinFrictionByAlphap_*alphaMax;
}

bool Foam::kineticTheoryModels::frictionalStressModels::Schaeffer::read()
{
    phi_.read(coeffDict());
    phi_ *= constant::mathematical::pi/180.0;
    alphaMinFrictionByAlphap_.read(coeffDict());

    return true;
}


// ************************************************************************* //
