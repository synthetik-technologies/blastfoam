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

#include "Wellek.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace aspectRatioModels
{
    defineTypeNameAndDebug(Wellek, 0);
    addToRunTimeSelectionTable
    (
        aspectRatioModel,
        Wellek,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::aspectRatioModels::Wellek::Wellek
(
    const dictionary& dict,
    const phasePair& pair
)
:
    aspectRatioModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::aspectRatioModels::Wellek::~Wellek()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::aspectRatioModels::Wellek::E
(
    const label nodei,
    const label nodej
) const
{
    return scalar(1)/(scalar(1) + 0.163*pow(pair_.Eo(nodei,nodej), 0.757));
}



Foam::scalar Foam::aspectRatioModels::Wellek::E
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    return
        scalar(1)/(scalar(1) + 0.163*pow(pair_.Eo(celli, nodei,nodej), 0.757));
}


// ************************************************************************* //
