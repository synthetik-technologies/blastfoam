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
    const kineticTheorySystem& kt
)
:
    radialModel(dict, kt),
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
    const volScalarField& alphap = kt_.alphap();
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

    forAll(kt_.phaseIndexes(), phaseI)
    {
        const phaseModel& phase =
            kt_.fluid().phases()[kt_.phaseIndexes()[phaseI]];
        alphard += volScalarField(phase)/phase.d();
    }

    return
        1.0/max(alphag, residualAlpha_)
      + 3.0*phase1.d()*phase2.d()*alphard
       /(sqr(max(alphag, residualAlpha_))*(phase1.d() + phase2.d()));
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::Lebowitz::gs0prime
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    const volScalarField& alphap = kt_.alphap();
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

    forAll(kt_.phaseIndexes(), phaseI)
    {
        const phaseModel& phase =
            kt_.fluid().phases()[kt_.phaseIndexes()[phaseI]];
        if (phase.name() != phase1.name())
        {
            alphard += volScalarField(phase)/phase.d();
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


// ************************************************************************* //
