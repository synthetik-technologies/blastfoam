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

#include "Chao.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"
#include "kineticTheorySystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(Chao, 0);
    addToRunTimeSelectionTable(dragModel, Chao, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::Chao::Chao
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    dragModel(dict, pair, registerObject),
    kineticTheorySystem_
    (
        pair_.phase1().mesh().lookupObject<kineticTheorySystem>
        (
            "kineticTheorySystem"
        )
    ),
    e_(kineticTheorySystem_.es(pair))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::Chao::~Chao()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::Chao::CdRe
(
    const label nodei,
    const label nodej
) const
{
    FatalErrorInFunction
        << "Not implemented."
        << "Drag coefficient not defined for the Chao model."
        << exit(FatalError);

    return pair_.phase1();
}


Foam::tmp<Foam::volScalarField> Foam::dragModels::Chao::K
(
    const label nodei,
    const label nodej
) const
{
    const phaseModel& phase1 = pair_.phase1();
    const phaseModel& phase2 = pair_.phase2();

    const volScalarField& rho1(phase1.rho());
    const volScalarField& rho2(phase2.rho());
    const volScalarField& Theta1 = pair_.phase1().Theta()();
    const volScalarField& Theta2 = pair_.phase2().Theta()();
    const scalar pi(Foam::constant::mathematical::pi);

    tmp<volScalarField> gij(kineticTheorySystem_.gs0(phase1, phase2));
    volScalarField dij(0.5*(phase1.d(nodei) + phase2.d(nodej)));
    volScalarField m0
    (
        constant::mathematical::pi/6.0
       *(
            phase1.rho()*pow3(phase1.d(nodei))
          + phase2.rho()*pow3(phase2.d(nodej))
        )
    );
    volScalarField magUr(pair_.magUr(nodei, nodej));

    return
        phase1*phase2
       *rho1*rho2/m0*sqr(dij)*(1.0 + e_)*gij
       *(
            sqrt(2.0*pi)*(sqrt(Theta1) + sqrt(Theta2))
          - sqrt(2.0)*pow(Theta1*Theta2, 0.25)
          + 0.5*pi*magUr - 1.135*sqrt(magUr)
           *(
                pow(Theta1, 0.25) + pow(Theta2, 0.25)
              - 0.8*pow(Theta1*Theta2, 0.125)
            )
        );
}


Foam::tmp<Foam::surfaceScalarField> Foam::dragModels::Chao::Kf
(
    const label nodei,
    const label nodej
) const
{
    return fvc::interpolate(K(nodei, nodej));
}


Foam::scalar Foam::dragModels::Chao::CdRei
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    FatalErrorInFunction
        << "Not implemented."
        << "Drag coefficient not defined for the Chao model."
        << exit(FatalError);

    return 0.0;
}


Foam::scalar Foam::dragModels::Chao::Ki
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    const phaseModel& phase1 = pair_.phase1();
    const phaseModel& phase2 = pair_.phase2();

    scalar rho1(phase1.rho()[celli]);
    scalar rho2(phase2.rho()[celli]);
    scalar Theta1(pair_.phase1().Theta()()[celli]);
    scalar Theta2(pair_.phase2().Theta()()[celli]);

    const scalar pi(Foam::constant::mathematical::pi);

    scalar gij(kineticTheorySystem_.gs0(phase1, phase2)()[celli]);
    scalar dij
    (
        0.5*(phase1.di(celli, nodei) + phase2.di(celli, nodej))
    );
    scalar m0
    (
        constant::mathematical::pi/6.0
       *(
            rho1*pow3(phase1.di(celli, nodei))
          + rho2*pow3(phase2.di(celli, nodej))
        )
    );
    scalar magUr(pair_.magUri(celli, nodei, nodej));

    return
        phase1[celli]*phase2[celli]
       *rho1*rho2/m0*sqr(dij)*(1.0 + e_)*gij
       *(
            sqrt(2.0*pi)*(sqrt(Theta1) + sqrt(Theta2))
          - sqrt(2.0)*pow(Theta1*Theta2, 0.25)
          + 0.5*pi*magUr - 1.135*sqrt(magUr)
           *(
                pow(Theta1, 0.25) + pow(Theta2, 0.25)
              - 0.8*pow(Theta1*Theta2, 0.125)
            )
        );
}

// ************************************************************************* //
