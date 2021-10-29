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

#include "movingWallImmersedBoundaryFvPatchVectorField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingWallImmersedBoundaryFvPatchVectorField::
movingWallImmersedBoundaryFvPatchVectorField
(
    volVectorField& f,
    const dictionary& dict,
    const immersedBoundaryObject& ibo
)
:
    immersedBoundaryVectorPatchField(f, dict, ibo)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::movingWallImmersedBoundaryFvPatchVectorField::updateCoeffs() const
{
    values_ = ibm_.velocity();
}


void Foam::movingWallImmersedBoundaryFvPatchVectorField::setValues()
{
    vector avg(gSum(values_*ibm_.magSf())/gSum(ibm_.magSf()));
    this->ibm_.setInternal
    (
        field_,
        avg//this->ibm_.velocity(ibm_.internalC())()
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeImmersedPatchTypeField
    (
        immersedBoundaryVectorPatchField,
        movingWallImmersedBoundaryFvPatchVectorField
    );
}
// ************************************************************************* //
