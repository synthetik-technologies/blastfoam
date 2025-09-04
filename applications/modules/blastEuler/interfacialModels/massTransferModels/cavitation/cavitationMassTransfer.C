/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derived work of OpenFOAM.

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

#include "cavitationMassTransfer.H"
#include "phaseSystem.H"
#include "phasePair.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace massTransferModels
{
    defineTypeNameAndDebug(cavitation, 0);
    addToRunTimeSelectionTable(massTransferModel, cavitation, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massTransferModels::cavitation::cavitation
(
    const dictionary& dict,
    const phasePair& pair
)
:
    massTransferModel(dict, pair),
    cavitation_
    (
        cavitationModel::New
        (
            dict,
            pair.phase1(),
            pair.phase1().fluidThermo(),
            pair.phase2(),
            pair.phase2().fluidThermo()
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::massTransferModels::cavitation::~cavitation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::massTransferModels::cavitation::K() const
{
    tmp<volScalarField> tResult =
        volScalarField::New
        (
            IOobject::groupName(typedName("dmdt"), pair_.name()),
            pair_.phase1().mesh(),
            dimDensity/dimTime,
            zeroGradientFvPatchField<scalar>::typeName
        );

    const Pair<tmp<volScalarField::Internal>> mDots(cavitation_->mDots());
    tResult.ref().internalFieldRef() = mDots[0] - mDots[1];
    tResult.ref().correctBoundaryConditions();

    return tResult;
}

// ************************************************************************* //
