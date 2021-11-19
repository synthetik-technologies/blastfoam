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

#include "symmetricTractionFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

symmetricTractionFvPatchVectorField::symmetricTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


symmetricTractionFvPatchVectorField::symmetricTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


symmetricTractionFvPatchVectorField::symmetricTractionFvPatchVectorField
(
    const symmetricTractionFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


symmetricTractionFvPatchVectorField::symmetricTractionFvPatchVectorField
(
    const symmetricTractionFvPatchVectorField& rifvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(rifvpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void symmetricTractionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const Field<vector>& rhoU =
        patch().lookupPatchField<surfaceVectorField, vector>("rhoUOwn");
    const Field<vector>& traction =
        patch().lookupPatchField<surfaceVectorField, vector>("tractionOwn");

    const Field<vector>& n =
        patch().lookupPatchField<surfaceVectorField, vector>("n");

    const Field<scalar>& pWaveSpeed =
        patch().lookupPatchField<volScalarField, scalar>("pWaveSpeed");

    Field<vector>::operator=
    (
        ((n*n) & (traction - pWaveSpeed*rhoU))
    );
    fixedValueFvPatchVectorField::updateCoeffs();
}


void symmetricTractionFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    symmetricTractionFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
