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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "solidWedgeFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidWedgeFvPatchVectorField::solidWedgeFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    wedgeFvPatchField<vector>(p, iF)
{}


solidWedgeFvPatchVectorField::solidWedgeFvPatchVectorField
(
    const solidWedgeFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    wedgeFvPatchField<vector>(ptf, p, iF, mapper)
{
    if (!isType<wedgeFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "solidWedgeFvPatchVectorField::solidWedgeFvPatchVectorField\n"
            "(\n"
            "    const solidWedgeFvPatchVectorField& ptf,\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<vector, volMesh>& iF,\n"
            "    const fvPatchFieldMapper& mapper\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }
}


solidWedgeFvPatchVectorField::solidWedgeFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    wedgeFvPatchField<vector>(p, iF, dict)
{
    if (!isType<wedgeFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "solidWedgeFvPatchVectorField::solidWedgeFvPatchVectorField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const Field<vector>& field,\n"
            "    dictionary& dict\n"
            ")\n",
            dict
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }

    //evaluate();

    // to initialise, we evaluate without non-orthogonal correction
    // this is because the explicit gradient is need for the correction
    // and obviously cannot be created until the field is created
    fvPatchField<vector>::operator==
    (
        transform
        (
            refCast<const wedgeFvPatch>(this->patch()).faceT(),
            this->patchInternalField()
        )
    );
}


solidWedgeFvPatchVectorField::solidWedgeFvPatchVectorField
(
    const solidWedgeFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    wedgeFvPatchField<vector>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<Field<vector> > solidWedgeFvPatchVectorField::snGrad() const
{
    //Info<< "solidWedgeFvPatchVectorField::snGrad()" << endl;
    // Method:
    // We project the patch normal back to the wedge centrePlane and find the
    // intersection point. We then extrapolate U from the Cn to this
    // intersection point (projC) using gradU looked up from the solver.
    // We can then calculate snGrad by getting the difference between projU and
    // transformed projU (i.e. across the wedge), and then dividing by the
    // magnitude of the distance between them.

    const wedgePolyPatch& wedgePatch =
        refCast<const wedgePolyPatch>(patch().patch());

    const vectorField& patchC = patch().patch().faceCentres();
    vectorField nHat(this->patch().nf());
    const vector& centreN = wedgePatch.centreNormal();
    scalarField d(((patch().Cn() - patchC) & centreN)/(nHat & centreN));
    vectorField projC(d*nHat + patchC);

    // Calculate correction vector which connects actual cell centre to the
    // transformed cell centre
    const vectorField k(projC - patch().Cn());

    Field<vector> pif(this->patchInternalField());

    const fvPatchField<tensor>& gradU =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + this->internalField().name() + ")"
        );

    Field<vector> projU
    (
        this->patchInternalField() + (k & gradU.patchInternalField())
    );

    // Calculate delta coeffs from proj position on centre plane to transformed
    // projected position
    scalarField projDeltaCoeff
    (
        1.0/mag(transform(wedgePatch.cellT(), projC) - projC)
    );

    return
    (
        transform(wedgePatch.cellT(), projU) - projU
    )*projDeltaCoeff;

    // old way without correction
    // return
    // (
    //  transform(refCast<const wedgeFvPatch>(this->patch()).cellT(), pif) - pif
    // )*(0.5*this->patch().deltaCoeffs());
}


void solidWedgeFvPatchVectorField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // Standard non-orthogonal correction vectors are parallel to the patch, but
    // for axisymmetric cases we must use a correction vector that lies along
    // the centre plane. Then we can transform/rotate this centre plane value to
    // the wedge patch.
    const wedgeFvPatch& wedgePatch =
        refCast<const wedgeFvPatch>(this->patch());

    // Rotate patchC field back to centre plane to find transformed cell centres
    const vectorField& patchC = patch().patch().faceCentres();
    vectorField transC(wedgePatch.faceT().T() & patchC);

    // Calculate correction vector which connects actual cell centre to the
    // transformed cell centre
    const vectorField k(transC - patch().Cn());

    const fvPatchField<tensor>& gradU =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + this->internalField().name() + ")"
        );

    Field<vector> pif(this->patchInternalField());
    pif += (k & gradU.patchInternalField());

    Field<vector>::operator=
    (
        transform(wedgePatch.faceT(), pif)
    );
}


void solidWedgeFvPatchVectorField::write(Ostream& os) const
{
    wedgeFvPatchVectorField::write(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, solidWedgeFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
