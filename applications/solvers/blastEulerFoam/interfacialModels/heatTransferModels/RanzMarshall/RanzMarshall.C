/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

#include "RanzMarshall.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(RanzMarshall, 0);
    addToRunTimeSelectionTable(heatTransferModel, RanzMarshall, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::RanzMarshall::RanzMarshall
(
    const dictionary& dict,
    const phasePair& pair
)
:
    heatTransferModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatTransferModels::RanzMarshall::~RanzMarshall()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::heatTransferModels::RanzMarshall::K
(
    const label nodei,
    const label nodej
) const
{
    const volScalarField& alphag(pair_.continuous().volumeFraction(nodej));
    volScalarField Pr(pair_.Pr(nodei, nodej));
    volScalarField Re(pair_.Re(nodei, nodej));
    volScalarField Nu
    (
        (7.0 - 10.0*alphag + 5.0*sqr(alphag))
       *(1.0 + 0.7*pow(Re, 0.2)*pow(Pr, 1.0/3.0))
      + (1.33 - 2.4*alphag + 1.2*sqr(alphag))
       *pow(Re, 0.7)*pow(Pr, 1.0/3.0)
//         scalar(2)
//       + 0.6*sqrt(pair_.Re(nodei, nodej))*cbrt(pair_.Pr(nodei, nodej))
    );

    return
        6.0
       *max(pair_.dispersed().volumeFraction(nodei), residualAlpha_)
       *pair_.continuous().kappa()
       *Nu
       /sqr(pair_.dispersed().d(nodei));
}


// ************************************************************************* //
