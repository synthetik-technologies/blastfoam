/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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

#include "PrincetonFrictionalStress.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace frictionalStressModels
{
    defineTypeNameAndDebug(Princeton, 0);

    addToRunTimeSelectionTable
    (
        frictionalStressModel,
        Princeton,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::Princeton::
Princeton
(
    const dictionary& dict,
    const kineticTheorySystem& kt
)
:
    frictionalStressModel(dict, kt),
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

Foam::kineticTheoryModels::frictionalStressModels::Princeton::
~Princeton()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::Princeton::
frictionalPressure
(
    const phaseModel& phase,
    const volScalarField& alphap,
    const volScalarField& alphaMax
) const
{
    volScalarField alphaMinSchaefer(alphaMinFrictionByAlphap_*alphaMax);

    return
        neg(alphap - alphaMinSchaefer)
       *Fr_*pow(max(alphap - alphaMinFriction_, scalar(0)), eta_)
       /pow(max(alphaMax - alphap, alphaDeltaMin_), p_)
      + dimensionedScalar("1e24", dimensionSet(1, -1, -2, 0, 0), 1e24)
       *pow(Foam::max(alphap - alphaMinSchaefer, scalar(0)), 10.0);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::Princeton::
frictionalPressurePrime
(
    const phaseModel& phase,
    const volScalarField& alphap,
    const volScalarField& alphaMax
) const
{
    volScalarField alphaMinSchaefer(alphaMinFrictionByAlphap_*alphaMax);

    return
        neg(alphap - alphaMinSchaefer)
       *Fr_
       *(
            eta_*pow(max(alphap - alphaMinFriction_, scalar(0)), eta_ - 1.0)
           *(alphaMax - alphap)
          + p_*pow(max(alphap - alphaMinFriction_, scalar(0)), eta_)
        )/pow(max(alphaMax - alphap, alphaDeltaMin_), p_ + 1.0)
      + dimensionedScalar("1e25", dimensionSet(1, -1, -2, 0, 0), 1e25)
       *pow(Foam::max(alphap - alphaMinSchaefer, scalar(0)), 9.0);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::Princeton::mu
(
    const phaseModel& phase,
    const volScalarField& alphap,
    const volScalarField& alphaMax,
    const volScalarField& Pc
) const
{
    volScalarField alphaMinFriction(alphaMinFrictionByAlphap_*alphaMax);
    tmp<volScalarField> da = phase.d();
    const volScalarField& Theta =
        phase.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName("Theta", alphap.group())
        );

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

    volScalarField divU(fvc::div(phase.phi()));
    volSymmTensorField S(dev(symm(fvc::grad(phase.U()))));
    volScalarField Sdd(S && S);
    volScalarField n
    (
        sqrt(3.0)/2.0*sin(phi_)*pos(divU)
      + 1.03*neg(divU)
    );

    tmp<volScalarField> PfByPc
    (
        pow
        (
            max
            (
                1.0
              - divU
               /max
                (
                    n*sqrt(2.0)*sin(phi_)*sqrt(Sdd + Theta/sqr(da())),
                    dimensionedScalar("small", divU.dimensions(), small)
                ),
                1e-6
            ),
            n - 1
        )
       *phase/max(alphap, phase.residualAlpha())
    );
    PfByPc.ref().max(1e-10);
    tmp<volScalarField> Pf(PfByPc()*Pc);

    forAll(Pc, celli)
    {
        if (alphap[celli] > alphaMinFriction[celli])
        {
            muf[celli] =
                sqrt(2.0)*Pf()[celli]*sin(phi_.value())
               /sqrt(Sdd[celli] + Theta[celli]/sqr(da()[celli]))
               *(
                    n[celli]
                  - (n[celli] - 1.0)
                   *pow(PfByPc()[celli], 1.0/(n[celli] - 1.0))
                );
        }
    }

    const fvPatchList& patches = phase.mesh().boundary();

    volScalarField::Boundary& mufBf = muf.boundaryFieldRef();

    forAll(patches, patchi)
    {
        if (!patches[patchi].coupled())
        {
            mufBf[patchi] =
            (
                0.5
               *Pc.boundaryField()[patchi]
               *sin(phi_.value())
            );
        }
    }

    // Correct coupled BCs
    muf.correctBoundaryConditions();

    return tmu;
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::Princeton::
alphaMinFriction
(
    const volScalarField& alphap,
    const volScalarField& alphaMax
) const
{
    volScalarField alphaMinSchaefer(alphaMinFrictionByAlphap_*alphaMax);

    return
        pos0(alphap - alphaMinSchaefer)*alphaMinSchaefer
      + neg(alphap - alphaMinSchaefer)*alphaMinFriction_;
}


bool Foam::kineticTheoryModels::frictionalStressModels::Princeton::read()
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
