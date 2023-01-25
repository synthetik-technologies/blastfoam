/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "LebowitzRadial.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace radialModels
{
    defineTypeNameAndDebug(Lebowitz, 0);

    addToRunTimeSelectionTable
    (
        radialModel,
        Lebowitz,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModels::Lebowitz::Lebowitz
(
    const dictionary& dict,
    const masterSystem& system
)
:
    radialModel(dict, system),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        dict.lookup("residualAlpha")
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModels::Lebowitz::~Lebowitz()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::Lebowitz::gs0
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    const volScalarField& alphap = system_.alpha();
    volScalarField alphag(1.0 - alphap);
    volScalarField alphard
    (
        IOobject
        (
            "alphard",
            alphap.time().timeName(),
            alphap.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        alphap.mesh(),
        dimensionedScalar("0", inv(dimLength), 0.0)
    );

    forAll(system_.phaseIndexes(), phaseI)
    {
        const phaseModel& phase =
            system_.fluid().phases()[system_.phaseIndexes()[phaseI]];
        alphard += volScalarField(phase)/phase.d();
    }

    return
        1.0/max(alphag, residualAlpha_)
      + 3.0*phase1.d()*phase2.d()*alphard
       /(sqr(max(alphag, residualAlpha_))*(phase1.d() + phase2.d()));
}


Foam::scalar
Foam::kineticTheoryModels::radialModels::Lebowitz::cellgs0
(
    const label celli,
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    const scalar alphap = system_.alpha()[celli];
    scalar alphag(1.0 - alphap);
    scalar alphard = 0.0;
    forAll(system_.phaseIndexes(), phaseI)
    {
        const phaseModel& phase =
            system_.fluid().phases()[system_.phaseIndexes()[phaseI]];
        alphard += phase[celli]/phase.celld(celli);
    }

    return
        1.0/max(alphag, residualAlpha_.value())
      + 3.0*phase1.celld(celli)*phase2.celld(celli)*alphard
       /(
           sqr(max(alphag, residualAlpha_.value()))
          *(phase1.celld(celli) + phase2.celld(celli))
        );
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::Lebowitz::gs0prime
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    const volScalarField& alphap = system_.alpha();
    volScalarField alphag(1.0 - alphap);
    volScalarField alphard
    (
        IOobject
        (
            "alphard",
            alphap.time().timeName(),
            alphap.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        alphap.mesh(),
        dimensionedScalar("0", inv(dimLength), 0.0)
    );

    forAll(system_.phaseIndexes(), phaseI)
    {
        const phaseModel& phase =
            system_.fluid().phases()[system_.phaseIndexes()[phaseI]];
        if (phase.name() != phase1.name())
        {
            alphard += phase/phase.d();
        }
    }

    volScalarField d1(phase1.d());
    volScalarField d2(phase2.d());
    return
        1.0/max(sqr(alphag), residualAlpha_)
       *(
            1.0
          + 3.0*d2/(d1 + d2)
           *(2.0/max(alphag, residualAlpha_)*(d1*alphard + phase1) + 1.0)
        );
}


Foam::scalar
Foam::kineticTheoryModels::radialModels::Lebowitz::cellgs0prime
(
    const label celli,
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    const scalar alphap = system_.alpha()[celli];
    scalar alphag(1.0 - alphap);
    scalar alphard = 0.0;

    forAll(system_.phaseIndexes(), phaseI)
    {
        const phaseModel& phase =
            system_.fluid().phases()[system_.phaseIndexes()[phaseI]];
        if (&phase != &phase1)
        {
            alphard += phase[celli]/phase.celld(celli);
        }
    }

    scalar d1(phase1.celld(celli));
    scalar d2(phase2.celld(celli));
    return
        1.0/max(sqr(alphag), residualAlpha_.value())
       *(
            1.0
          + 3.0*d2/(d1 + d2)
           *(
               2.0/max(alphag, residualAlpha_.value())
              *(d1*alphard + phase1[celli])
             + 1.0
            )
        );
}


// ************************************************************************* //
