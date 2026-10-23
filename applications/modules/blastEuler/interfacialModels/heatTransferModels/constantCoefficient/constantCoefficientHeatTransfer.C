/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2025-06-09 Jeff Heylmun:    Added cell based returns
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

#include "constantCoefficientHeatTransfer.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(constantCoefficient, 0);
    addToRunTimeSelectionTable(heatTransferModel, constantCoefficient, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::constantCoefficient::constantCoefficient
(
    const dictionary& dict,
    const phasePair& pair
)
:
    dispersedHeatTransferModel(dict, pair),
    Kh_("Kh", dimK, dict.lookup<scalar>("Kh")),
    fieldName_(dict.lookupOrDefault<word>("field", word::null))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatTransferModels::constantCoefficient::~constantCoefficient()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volScalarField*
Foam::heatTransferModels::constantCoefficient::field() const
{
    if (fieldName_ == word::null)
    {
        return nullptr;
    }
    return &pair_.phase1().mesh().lookupObject<volScalarField>(fieldName_);
}


Foam::tmp<Foam::volScalarField>
Foam::heatTransferModels::constantCoefficient::Ki() const
{
    const volScalarField* fPtr = field();
    if (fPtr)
    {
        return volScalarField::New(typeName + ":Ki", Kh_*max(*fPtr, scalar(0)));
    }
    return volScalarField::New
    (
        typeName + ":Ki",
        pair_.phase1().mesh(),
        Kh_
    );
}


Foam::scalar Foam::heatTransferModels::constantCoefficient::cellKi
(
    const label celli
) const
{
    const volScalarField* fPtr = field();
    if (fPtr)
    {
        return Kh_.value()*max((*fPtr)[celli], scalar(0));
    }
    return Kh_.value();
}


// ************************************************************************* //
