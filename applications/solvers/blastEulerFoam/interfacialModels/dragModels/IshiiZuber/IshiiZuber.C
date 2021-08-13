/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-18 Jeff Heylmun:    Added support of polydisperse phase models
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

#include "IshiiZuber.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(IshiiZuber, 0);
    addToRunTimeSelectionTable(dragModel, IshiiZuber, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::IshiiZuber::IshiiZuber
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    dragModel(dict, pair, registerObject)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::IshiiZuber::~IshiiZuber()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::dragModels::IshiiZuber::CdRe
(
    const label nodei,
    const label nodej
) const
{
    volScalarField Re(pair_.Re(nodei, nodej));
    volScalarField Eo(pair_.Eo(nodei, nodej));

    volScalarField mud(pair_.dispersed().mu());
    volScalarField muc(pair_.continuous().mu());

    volScalarField muStar((mud + 0.4*muc)/(mud + muc));

    volScalarField muMix
    (
        muc
       *pow
        (
            max(1 - pair_.dispersed().volumeFraction(nodei),
            scalar(1e-3)), -2.5*muStar
        )
    );

    volScalarField ReM(Re*muc/muMix);
    volScalarField CdRe
    (
        pos0(1000 - ReM)*24.0*(scalar(1) + 0.15*pow(ReM, 0.687))
      + neg(1000 - ReM)*0.44*ReM
    );

    volScalarField F
    (
        (muc/muMix)*sqrt(1 - pair_.dispersed().volumeFraction(nodei))
    );
    F.max(1e-3);

    volScalarField Ealpha((1 + 17.67*pow(F, 0.8571428))/(18.67*F));

    volScalarField CdReEllipse(Ealpha*0.6666*sqrt(Eo)*Re);

    return
        pos0(CdReEllipse - CdRe)
       *min
        (
            CdReEllipse,
            Re*sqr(1 - pair_.dispersed().volumeFraction(nodej))*2.66667
        )
      + neg(CdReEllipse - CdRe)*CdRe;
}


Foam::scalar Foam::dragModels::IshiiZuber::CdRei
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    scalar Re(pair_.Rei(celli, nodei, nodej));
    scalar Eo(pair_.Eoi(celli, nodei, nodej));

    scalar mud(pair_.dispersed().nui(celli)*pair_.dispersed().rho()[celli]);
    scalar muc(pair_.continuous().nui(celli)*pair_.continuous().rho()[celli]);

    scalar muStar((mud + 0.4*muc)/(mud + muc));

    scalar muMix
    (
        muc
       *pow
        (
            max(1 - pair_.dispersed().volumeFractioni(celli, nodei),
            scalar(1e-3)), -2.5*muStar
        )
    );

    scalar ReM(Re*muc/muMix);
    scalar CdRe
    (
        pos0(1000 - ReM)*24.0*(scalar(1) + 0.15*pow(ReM, 0.687))
      + neg(1000 - ReM)*0.44*ReM
    );

    scalar F
    (
        (muc/muMix)*sqrt(1 - pair_.dispersed().volumeFractioni(celli, nodei))
    );
    F = max(F, 1e-3);

    scalar Ealpha((1 + 17.67*pow(F, 0.8571428))/(18.67*F));

    scalar CdReEllipse(Ealpha*0.6666*sqrt(Eo)*Re);

    return
        pos0(CdReEllipse - CdRe)
       *min
        (
            CdReEllipse,
            Re*sqr(1 - pair_.dispersed().volumeFractioni(celli, nodej))*2.66667
        )
      + neg(CdReEllipse - CdRe)*CdRe;
}

// ************************************************************************* //
