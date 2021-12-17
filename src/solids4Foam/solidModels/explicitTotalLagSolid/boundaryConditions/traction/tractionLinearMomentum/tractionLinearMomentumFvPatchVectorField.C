/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "tractionLinearMomentumFvPatchVectorField.H"
#include "solidTractionFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tractionLinearMomentumFvPatchVectorField::
tractionLinearMomentumFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


tractionLinearMomentumFvPatchVectorField::
tractionLinearMomentumFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict)
{}


tractionLinearMomentumFvPatchVectorField::
tractionLinearMomentumFvPatchVectorField
(
    const tractionLinearMomentumFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


tractionLinearMomentumFvPatchVectorField::
tractionLinearMomentumFvPatchVectorField
(
    const tractionLinearMomentumFvPatchVectorField& rifvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(rifvpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void tractionLinearMomentumFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const solidTractionFvPatchVectorField& pD =
        dynamicCast<const solidTractionFvPatchVectorField>
        (
            patch().lookupPatchField
            <
                volVectorField,
                vector
            >("D")
        );

    const fvsPatchField<vector>& n =
        patch().lookupPatchField<surfaceVectorField, vector>("n");

    const fvsPatchField<vector>& rhoU =
        patch().lookupPatchField<surfaceVectorField, vector>("rhoUOwn");

    const fvsPatchField<vector>& traction =
        patch().lookupPatchField<surfaceVectorField, vector>("tractionOwn");

    const Field<tensor>& stabTraction =
        patch().lookupPatchField<surfaceTensorField, tensor>("stabTraction");

    Field<vector>::operator=
    (
        rhoU
      + (
            stabTraction
          & (
                (pD.traction() - pD.pressure()*n) - traction
            )
        )
    );
    fixedValueFvPatchVectorField::updateCoeffs();
}


void tractionLinearMomentumFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    tractionLinearMomentumFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
