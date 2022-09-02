/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "hydrostaticPressureFvPatchVectorField.H"
#include "volFields.H"
#include "lookupSolidModel.H"
#include "uniformDimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

hydrostaticPressureFvPatchVectorField::
hydrostaticPressureFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(p, iF),
    pRef_(0.0),
    hRef_(0.0),
    rho_(1.0)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
    traction() = vector::zero;
    pressure() = 0.0;
}


hydrostaticPressureFvPatchVectorField::
hydrostaticPressureFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidTractionFvPatchVectorField(p, iF),
    pRef_(dict.lookup<scalar>("pRef")),
    hRef_(dict.lookup<scalar>("hRef")),
    rho_(dict.lookup<scalar>("rho"))
{
    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        fvPatchVectorField::operator=(patchInternalField());
    }

    traction() = vector::zero;
    pressure() = 0.0;
}


hydrostaticPressureFvPatchVectorField::
hydrostaticPressureFvPatchVectorField
(
    const hydrostaticPressureFvPatchVectorField& hpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidTractionFvPatchVectorField(hpvf, p, iF, mapper),
    pRef_(hpvf.pRef_),
    hRef_(hpvf.hRef_),
    rho_(hpvf.rho_)
{
    traction() = vector::zero;
    pressure() = 0.0;
}


hydrostaticPressureFvPatchVectorField::
hydrostaticPressureFvPatchVectorField
(
    const hydrostaticPressureFvPatchVectorField& hpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(hpvf, iF),
    pRef_(hpvf.pRef_),
    hRef_(hpvf.hRef_),
    rho_(hpvf.rho_)
{
    traction() = vector::zero;
    pressure() = 0.0;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void hydrostaticPressureFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const uniformDimensionedVectorField& g =
        this->patch().boundaryMesh().mesh().lookupObject
        <
            uniformDimensionedVectorField
        >("g");
    vector dir(g.value()/mag(g.value()));

    vectorField x(this->patch().Cf());
    if (!lookupSolidModel(this->patch().boundaryMesh().mesh()).movingMesh())
    {
        x += this->patch().lookupPatchField<volVectorField, vector>("D");
        if (internalField().name() == "DD")
        {
            x += this->patch().lookupPatchField<volVectorField, vector>("DD");
        }
    }
    else
    {
        x += this->patch().lookupPatchField<volVectorField, vector>("DD");
    }

    scalarField gh((x + dir*hRef_) & g.value());
    this->pressure() = pRef_ + rho_*gh;

    solidTractionFvPatchVectorField::updateCoeffs();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, hydrostaticPressureFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
