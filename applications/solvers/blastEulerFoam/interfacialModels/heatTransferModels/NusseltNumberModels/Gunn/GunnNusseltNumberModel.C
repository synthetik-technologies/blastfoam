/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technologies
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "GunnNusseltNumberModel.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace NusseltNumberModels
{
    defineTypeNameAndDebug(Gunn, 0);
    addToRunTimeSelectionTable(NusseltNumberModel, Gunn, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::NusseltNumberModels::Gunn::Gunn

(
    const dictionary& dict,
    const phasePair& pair
)
:
    NusseltNumberModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::NusseltNumberModels::Gunn::~Gunn()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::NusseltNumberModels::Gunn::Nu
(
    const label nodei,
    const label nodej
) const
{
    const volScalarField& alphag(pair_.continuous().volumeFraction(nodej));
    volScalarField Pr(pair_.Pr(nodei, nodej));
    volScalarField Re(pair_.Re(nodei, nodej));
    return
        (7.0 - 10.0*alphag + 5.0*sqr(alphag))
       *(1.0 + 0.7*pow(Re, 0.2)*pow(Pr, 1.0/3.0))
      + (1.33 - 2.4*alphag + 1.2*sqr(alphag))
       *pow(Re, 0.7)*pow(Pr, 1.0/3.0);
}


Foam::scalar Foam::NusseltNumberModels::Gunn::cellNu
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    scalar alphag(pair_.continuous().cellvolumeFraction(celli, nodej));
    scalar Pr(pair_.cellPr(celli, nodei, nodej));
    scalar Re(pair_.cellRe(celli, nodei, nodej));
    return
        (7.0 - 10.0*alphag + 5.0*sqr(alphag))
       *(1.0 + 0.7*pow(Re, 0.2)*pow(Pr, 1.0/3.0))
      + (1.33 - 2.4*alphag + 1.2*sqr(alphag))
       *pow(Re, 0.7)*pow(Pr, 1.0/3.0);
}


// ************************************************************************* //
