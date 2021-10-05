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

#include "ductNusseltNumberModel.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace NusseltNumberModels
{
    defineTypeNameAndDebug(duct, 0);
    addToRunTimeSelectionTable(NusseltNumberModel, duct, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::NusseltNumberModels::duct::duct
(
    const dictionary& dict,
    const phasePair& pair
)
:
    NusseltNumberModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::NusseltNumberModels::duct::~duct()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::NusseltNumberModels::duct::Nu
(
    const label nodei,
    const label nodej
) const
{
    tmp<volScalarField> Pr(pair_.Pr(nodei, nodej));
    tmp<volScalarField> Re(pair_.Re(nodei, nodej));
    return 0.023*pow(Re, 0.8)*pow(Pr, 0.4);
}


Foam::scalar Foam::NusseltNumberModels::duct::cellNu
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    scalar Pr(pair_.cellPr(celli, nodei, nodej));
    scalar Re(pair_.cellRe(celli, nodei, nodej));
    return 0.23*pow(Re, 0.8)*pow(Pr, 0.4);
}

// ************************************************************************* //
