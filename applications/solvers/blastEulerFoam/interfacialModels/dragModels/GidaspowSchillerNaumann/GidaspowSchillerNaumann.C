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

#include "GidaspowSchillerNaumann.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(GidaspowSchillerNaumann, 0);
    addToRunTimeSelectionTable(dragModel, GidaspowSchillerNaumann, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::GidaspowSchillerNaumann::GidaspowSchillerNaumann
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

Foam::dragModels::GidaspowSchillerNaumann::~GidaspowSchillerNaumann()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::dragModels::GidaspowSchillerNaumann::CdRe
(
    const label nodei,
    const label nodej
) const
{
    volScalarField alpha2
    (
        max
        (
            scalar(1) - pair_.dispersed().volumeFraction(),
            pair_.continuous().residualAlpha()
        )
    );

    volScalarField Re(alpha2*pair_.Re(nodei, nodej));

    volScalarField CdsRe
    (
        neg(Re - 1000)*24.0*(1.0 + 0.15*pow(Re, 0.687))/alpha2
      + pos0(Re - 1000)*0.44*max(Re, residualRe_)
    );

    return
        CdsRe
       *pow(alpha2, -2.65)
       *max
        (
            pair_.continuous().volumeFraction(nodej),
            pair_.continuous().residualAlpha()
        );
}


Foam::scalar Foam::dragModels::GidaspowSchillerNaumann::CdRei
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
            pair_.continuous()[celli],
            pair_.continuous().residualAlpha().value()
        )
    );

    scalar Re(alpha2*pair_.Rei(celli, nodei, nodej));

    scalar CdsRe
    (
        neg(Re - 1000)*24.0*(1.0 + 0.15*pow(Re, 0.687))/alpha2
      + pos0(Re - 1000)*0.44*max(Re, residualRe_.value())
    );

    return
        CdsRe
       *pow(alpha2, -2.65)
       *max
        (
            pair_.continuous()[celli],
            pair_.continuous().residualAlpha().value()
        );
}

// ************************************************************************* //
