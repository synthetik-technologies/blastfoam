/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "totalThermoPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::totalThermoPressureFvPatchScalarField::
totalThermoPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    dynamicThermoPressureFvPatchScalarField(p, iF),
    UName_(IOobject::groupName("U", this->phaseName_)),
    phiName_(IOobject::groupName("phi", this->phaseName_))
{}


Foam::totalThermoPressureFvPatchScalarField::
totalThermoPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    dynamicThermoPressureFvPatchScalarField(p, iF, dict),
    UName_
    (
        dict.lookupOrDefault<word>
        (
            "U",
            IOobject::groupName("U", this->phaseName_)
        )
    ),
    phiName_
    (
        dict.lookupOrDefault<word>
        (
            "phi",
            IOobject::groupName("phi", this->phaseName_)
        )
    )
{}


Foam::totalThermoPressureFvPatchScalarField::
totalThermoPressureFvPatchScalarField
(
    const totalThermoPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    dynamicThermoPressureFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_)
{}


Foam::totalThermoPressureFvPatchScalarField::
totalThermoPressureFvPatchScalarField
(
    const totalThermoPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    dynamicThermoPressureFvPatchScalarField(tppsf, iF),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::totalThermoPressureFvPatchScalarField::updateCoeffs
(
    const scalarField& p0p,
    const vectorField& Up
)
{
    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    dynamicThermoPressureFvPatchScalarField::updateCoeffs
    (
        p0_,
        0.5*(1 - pos0(phip))*magSqr(Up)
    );
}


void Foam::totalThermoPressureFvPatchScalarField::updateCoeffs()
{
    const fvPatchField<vector>& Up =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    updateCoeffs(p0_, Up);
}


void Foam::totalThermoPressureFvPatchScalarField::write(Ostream& os) const
{
    dynamicThermoPressureFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        totalThermoPressureFvPatchScalarField
    );
}

// ************************************************************************* //
