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

#include "GidaspowErgunWenYu.H"
#include "phasePair.H"
#include "Ergun.H"
#include "WenYu.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(GidaspowErgunWenYu, 0);
    addToRunTimeSelectionTable(dragModel, GidaspowErgunWenYu, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::GidaspowErgunWenYu::GidaspowErgunWenYu
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    dragModel(dict, pair, registerObject),
    Ergun_
    (
        new Ergun
        (
            dict,
            pair,
            false
        )
    ),
    WenYu_
    (
        new WenYu
        (
            dict,
            pair,
            false
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::GidaspowErgunWenYu::~GidaspowErgunWenYu()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::dragModels::GidaspowErgunWenYu::CdRe
(
    const label nodei,
    const label nodej
) const
{
    return
        pos0(pair_.continuous() - 0.8)*WenYu_->CdRe(nodei, nodej)
      + neg(pair_.continuous() - 0.8)*Ergun_->CdRe(nodei, nodej);
}


Foam::scalar Foam::dragModels::GidaspowErgunWenYu::CdRei
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    return
        pos0(pair_.continuous()[celli] - 0.8)
       *WenYu_->CdRei(celli, nodei, nodej)
      + neg(pair_.continuous()[celli] - 0.8)
       *Ergun_->CdRei(celli, nodei, nodej);
}

// ************************************************************************* //
