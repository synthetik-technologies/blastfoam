/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "fixedFaceFluxPressureFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedFaceFluxPressureFvPatchScalarField::fixedFaceFluxPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(p, iF),
    refFace_(0),
    refP_(0.0)
{}


Foam::fixedFaceFluxPressureFvPatchScalarField::fixedFaceFluxPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchScalarField(p, iF),
    refFace_(dict.lookupOrDefault("refFace", 0)),
    refP_(dict.lookupType<scalar>("refP"))
{
    fvPatchScalarField::operator=(this->patchInternalField());
}


Foam::fixedFaceFluxPressureFvPatchScalarField::fixedFaceFluxPressureFvPatchScalarField
(
    const fixedFaceFluxPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchScalarField(p, iF),
    refFace_(ptf.refFace_),
    refP_(ptf.refP_)
{}


Foam::fixedFaceFluxPressureFvPatchScalarField::fixedFaceFluxPressureFvPatchScalarField
(
    const fixedFaceFluxPressureFvPatchScalarField& wbppsf
)
:
    zeroGradientFvPatchScalarField(wbppsf),
    refFace_(wbppsf.refFace_),
    refP_(wbppsf.refP_)
{}


Foam::fixedFaceFluxPressureFvPatchScalarField::fixedFaceFluxPressureFvPatchScalarField
(
    const fixedFaceFluxPressureFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(wbppsf, iF),
    refFace_(wbppsf.refFace_),
    refP_(wbppsf.refP_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedFaceFluxPressureFvPatchScalarField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    scalarField internal(this->patchInternalField());
    if (Pstream::myProcNo() == Pstream::masterNo())
    {
        internal[refFace_] = refP_;
    }
    scalarField::operator=(internal);
    fvPatchScalarField::evaluate();
}

Foam::tmp<Foam::scalarField>
Foam::fixedFaceFluxPressureFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    tmp<scalarField> viCoeffs
    (
        new scalarField(this->size(), 1.0)
    );
    viCoeffs.ref()[refFace_] = 0.0;
    return viCoeffs;
}


Foam::tmp<Foam::scalarField>
Foam::fixedFaceFluxPressureFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    tmp<scalarField> vbCoeffs
    (
        new scalarField(this->size(), 0.0)
    );
    vbCoeffs.ref()[refFace_] = refP_;
    return vbCoeffs;
}

Foam::tmp<Foam::scalarField>
Foam::fixedFaceFluxPressureFvPatchScalarField::gradientInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    tmp<scalarField> viCoeffs
    (
        new scalarField(this->size(), 0.0)
    );
    viCoeffs.ref()[refFace_] = this->patch().deltaCoeffs()[refFace_];
    return viCoeffs;
}


Foam::tmp<Foam::scalarField>
Foam::fixedFaceFluxPressureFvPatchScalarField::gradientBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    tmp<scalarField> vbCoeffs
    (
        new scalarField(this->size(), 0.0)
    );
    vbCoeffs.ref()[refFace_] = this->patch().deltaCoeffs()[refFace_]*refP_;
    return vbCoeffs;
}

void Foam::fixedFaceFluxPressureFvPatchScalarField::write(Ostream& os) const
{
    zeroGradientFvPatchScalarField::write(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedFaceFluxPressureFvPatchScalarField
    );
}


// ************************************************************************* //
