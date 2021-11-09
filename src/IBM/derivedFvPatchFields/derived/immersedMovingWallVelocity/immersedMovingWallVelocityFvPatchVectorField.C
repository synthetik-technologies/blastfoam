/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technology
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

#include "immersedMovingWallVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedMovingWallVelocityFvPatchVectorField::immersedMovingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    immersedValueFvPatchField<vector>(p, iF)
{}


Foam::immersedMovingWallVelocityFvPatchVectorField::immersedMovingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    immersedValueFvPatchField<vector>(p, iF, dict, false)
{}


Foam::immersedMovingWallVelocityFvPatchVectorField::immersedMovingWallVelocityFvPatchVectorField
(
    const immersedMovingWallVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    immersedValueFvPatchField<vector>(ptf, p, iF, mapper)
{}


Foam::immersedMovingWallVelocityFvPatchVectorField::immersedMovingWallVelocityFvPatchVectorField
(
    const immersedMovingWallVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    immersedValueFvPatchField<vector>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::immersedMovingWallVelocityFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    this->immersedField_ = this->object_.velocity();
    immersedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::immersedMovingWallVelocityFvPatchVectorField::setInternal()
{
    if (!this->setInternal_)
    {
        return;
    }
    const labelList internalCells(this->object_.internalCells());
    vectorField points
    (
        UIndirectList<vector>
        (
            this->internalField().mesh().C(),
            internalCells
        )()
    );
    vectorField v(this->object_.velocity(points));
    this->object_.setInternal(this->internalFieldRef(), v);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        immersedMovingWallVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
