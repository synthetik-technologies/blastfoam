/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "WenYu.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(WenYu, 0);
    addToRunTimeSelectionTable(dragModel, WenYu, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::WenYu::WenYu
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    dragModel(dict, pair, registerObject),
    residualRe_("residualRe", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::WenYu::~WenYu()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::WenYu::CdRe
(
    const label nodei,
    const label nodej
) const
{
    volScalarField alpha2
    (
        max
        (
            1.0 - pair_.dispersed().volumeFraction(nodei),
            pair_.continuous().residualAlpha()
        )
    );

    volScalarField Res(alpha2*pair_.Re(nodei, nodej));
    volScalarField CdsRes
    (
        neg(Res - 1000)*24.0*(1.0 + 0.15*pow(Res, 0.687))
      + pos0(Res - 1000)*0.44*max(Res, residualRe_)
    );

    return
        CdsRes
       *pow(alpha2, -3.65)
       *max(pair_.continuous(), pair_.continuous().residualAlpha());
}


Foam::scalar Foam::dragModels::WenYu::CdRe
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    scalar alpha2
    (
        max
        (
            1.0 - pair_.dispersed().volumeFractioni(celli, nodei),
            pair_.continuous().residualAlpha().value()
        )
    );

    scalar Res(alpha2*pair_.Re(celli, nodei, nodej));
    scalar CdsRes
    (
        neg(Res - 1000)*24.0*(1.0 + 0.15*pow(Res, 0.687))
      + pos0(Res - 1000)*0.44*max(Res, residualRe_.value())
    );

    return
        CdsRes
       *pow(alpha2, -3.65)
       *max
        (
            pair_.continuous()[celli],
            pair_.continuous().residualAlpha().value()
        );
}

// ************************************************************************* //
