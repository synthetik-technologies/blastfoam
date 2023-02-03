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
    rho_(1.0),
    phFunc_(nullptr)
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
    pRef_(0.0),
    hRef_(0.0),
    rho_(1.0),
    phFunc_(nullptr)
{
    if (dict.found("phFunc"))
    {
        phFunc_ = Function1<scalar>::New("phFunc", dict);
    }
    else
    {
        dict.lookup("pRef") >> pRef_;
        dict.lookup("hRef") >> hRef_;
        dict.lookup("rho") >> rho_;
    }

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
    rho_(hpvf.rho_),
    phFunc_(hpvf.phFunc_, false)
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
    rho_(hpvf.rho_),
    phFunc_(hpvf.phFunc_, false)
{
    traction() = vector::zero;
    pressure() = 0.0;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool hydrostaticPressureFvPatchVectorField::updateFields()
{
    const uniformDimensionedVectorField& g =
        this->patch().boundaryMesh().mesh().lookupObject
        <
            uniformDimensionedVectorField
        >("g");

    vectorField x(this->patch().Cf());
    if (!lookupSolidModel(this->patch().boundaryMesh().mesh()).movingMesh())
    {
        x += this->patch().lookupPatchField<volVectorField, vector>("D");
    }
    else
    {
        x += this->patch().lookupPatchField<volVectorField, vector>("DD");
    }

    if (phFunc_.valid())
    {
        scalarField h((x & g.value())/mag(g.value()) + hRef_);
        forAll(h, i)
        {
            h[i] = phFunc_->value(h[i]);
        }
        this->pressure() = h;
    }
    else
    {
        scalarField gh((x & g.value()) + mag(g.value())*hRef_);
        this->pressure() = pRef_ + rho_*gh;
    }

    return false;
}


void Foam::hydrostaticPressureFvPatchVectorField::write(Ostream& os) const
{
    solidTractionFvPatchVectorField::write(os);
    writeEntry(os, "pRef", pRef_);
    writeEntry(os, "hRef", hRef_);
    writeEntry(os, "rho", rho_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, hydrostaticPressureFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
