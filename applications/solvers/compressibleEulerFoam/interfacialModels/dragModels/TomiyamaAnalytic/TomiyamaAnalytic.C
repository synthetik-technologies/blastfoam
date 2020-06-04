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

#include "TomiyamaAnalytic.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(TomiyamaAnalytic, 0);
    addToRunTimeSelectionTable(dragModel, TomiyamaAnalytic, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::TomiyamaAnalytic::TomiyamaAnalytic
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    dragModel(dict, pair, registerObject),
    residualRe_("residualRe", dimless, dict),
    residualEo_("residualEo", dimless, dict),
    residualE_("residualE", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::TomiyamaAnalytic::~TomiyamaAnalytic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::dragModels::TomiyamaAnalytic::CdRe
(
    const label nodei,
    const label nodej
) const
{
    volScalarField Eo(max(pair_.Eo(nodei, nodej), residualEo_));
    volScalarField E(max(pair_.E(nodei, nodej), residualE_));

    volScalarField OmEsq(max(scalar(1) - sqr(E), sqr(residualE_)));
    volScalarField rtOmEsq(sqrt(OmEsq));

    volScalarField F(max(asin(rtOmEsq) - E*rtOmEsq, residualE_)/OmEsq);

    return
        (8.0/3.0)
       *Eo
       /(
            Eo*pow(E, 2.0/3.0)/OmEsq
          + 16*pow(E, 4.0/3.0)
        )
       /sqr(F)
       *max(pair_.Re(nodei, nodej), residualRe_);
}


Foam::scalar Foam::dragModels::TomiyamaAnalytic::CdRe
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    scalar Eo(max(pair_.Eo(celli, nodei, nodej), residualEo_.value()));
    scalar E(max(pair_.E(celli, nodei, nodej), residualE_.value()));

    scalar OmEsq(max(scalar(1) - sqr(E), sqr(residualE_.value())));
    scalar rtOmEsq(sqrt(OmEsq));

    scalar F(max(asin(rtOmEsq) - E*rtOmEsq, residualE_.value())/OmEsq);

    return
        (8.0/3.0)
       *Eo
       /(
            Eo*pow(E, 2.0/3.0)/OmEsq
          + 16*pow(E, 4.0/3.0)
        )
       /sqr(F)
       *max(pair_.Re(celli, nodei, nodej), residualRe_.value());
}

// ************************************************************************* //
