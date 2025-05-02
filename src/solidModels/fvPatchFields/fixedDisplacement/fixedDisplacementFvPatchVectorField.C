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

#include "fixedDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "valuePointPatchFields.H"
#include "movingObject.H"
#include "uniformDimensionedFields.H"


// * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * * //

void Foam::fixedDisplacementFvPatchVectorField::makeInterp() const
{
    if (interpPtr_.valid())
    {
        FatalErrorInFunction
            << "pointer already set" << endl
            << abort(FatalError);
    }

    interpPtr_.set(new primitivePatchInterpolation(patch().patch()));
}


const Foam::fixedDisplacementFvPatchVectorField::primitivePatchInterpolation&
Foam::fixedDisplacementFvPatchVectorField::interp() const
{
    if (interpPtr_.empty())
    {
        makeInterp();
    }

    return interpPtr_();
}


void Foam::fixedDisplacementFvPatchVectorField::setPointDisplacement
(
    const vectorField& faceDisp
)
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if (mesh.foundObject<pointVectorField>("point" + internalField().name()))
    {
        const pointVectorField& pointD =
            mesh.lookupObject<pointVectorField>
            (
                "point" + internalField().name()
            );

        // Check if the boundary is fixedValue
        if
        (
            isA<valuePointPatchVectorField>
            (
                pointD.boundaryField()[patch().index()]
            )
        )
        {
            // Use const_cast to set boundary condition
            valuePointPatchVectorField& patchPointD =
                refCast<valuePointPatchVectorField>
                (
                    const_cast<pointVectorField&>
                    (
                        pointD
                    ).boundaryFieldRef()[patch().index()]
                );

            // Interpolate face values to the points
            patchPointD == interp().faceToPointInterpolate(faceDisp);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    totalDisp_(p.size(), vector::zero),
    dispSeries_(),
    displacementName_(word::null),
    regionName_(word::null),
    interpPtr_()
{}


Foam::fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fixedDisplacementFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    totalDisp_(mapper(ptf.totalDisp_)),
    dispSeries_(ptf.dispSeries_, false),
    displacementName_(word::null),
    regionName_(word::null),
    interpPtr_()
{}


Foam::fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    totalDisp_(*this),
    dispSeries_(),
    displacementName_(word::null),
    regionName_(word::null),
    interpPtr_()
{
    // Check if displacement is time-varying
    if (dict.found("displacementSeries"))
    {
        DebugInfo
            << "    displacement is time-varying" << endl;
        dispSeries_ =
            Function1<vector>::New
            (
                "displacementSeries",
                this->db().time().userUnits(),
                dimLength,
                dict
            );

        Field<vector>::operator=
        (
            dispSeries_->value(this->db().time().value())
        );
    }
    else if (dict.found("displacementName"))
    {
        dict.readIfPresent("displacementName", displacementName_);
        dict.readIfPresent("displacementRegion", regionName_);
    }

}


Foam::fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fixedDisplacementFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    totalDisp_(pivpvf.totalDisp_),
    dispSeries_(pivpvf.dispSeries_, false),
    displacementName_(pivpvf.displacementName_),
    regionName_(pivpvf.regionName_),
    interpPtr_()
{}


// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //


Foam::fixedDisplacementFvPatchVectorField::~fixedDisplacementFvPatchVectorField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedDisplacementFvPatchVectorField::map
(
    const fvPatchField<vector>& ptf,
    const fieldMapper& mapper
)
{
    fixedValueFvPatchVectorField::map(ptf, mapper);

    const fixedDisplacementFvPatchVectorField& fdpvf =
        refCast<const fixedDisplacementFvPatchVectorField>(ptf);

    mapper(totalDisp_, fdpvf.totalDisp_);
    interpPtr_.clear();
}


void Foam::fixedDisplacementFvPatchVectorField::reset
(
    const fvPatchField<vector>& ptf
)
{
    fixedValueFvPatchVectorField::reset(ptf);

    const fixedDisplacementFvPatchVectorField& fdpvf =
        refCast<const fixedDisplacementFvPatchVectorField>(ptf);

    totalDisp_.reset(fdpvf.totalDisp_);
    interpPtr_.clear();
}


void Foam::fixedDisplacementFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (dispSeries_.valid())
    {
        vectorField disp
        (
            this->size(),
            dispSeries_->value(this->db().time().value())
        );
        if (internalField().name() == "DD")
        {
            // Incremental approach, so we wil set the increment of displacement
            // Lookup the old displacement field and subtract it from the total
            // displacement
            const volVectorField& Dold =
                db().lookupObject<volVectorField>("D").oldTime();

            disp -= Dold.boundaryField()[patch().index()];
        }
        Field<vector>::operator=(disp);
    }
    else if (!displacementName_.empty())
    {
        const objectRegistry& obr =
            !regionName_.empty()
          ? db().time().lookupObject<objectRegistry>(regionName_)
          : db();
        if (obr.foundObject<movingObject>(displacementName_))
        {
            const movingObject& object =
                obr.lookupObject<movingObject>(displacementName_);
            Field<vector>::operator=
            (
                object.centreOfMass() - object.initialCentreOfMass()
            );
        }
        else if
        (
            obr.foundObject<uniformDimensionedVectorField>(displacementName_)
        )
        {
            const uniformDimensionedVectorField& disp =
                obr.lookupObject<uniformDimensionedVectorField>
                (
                    displacementName_
                );
            Field<vector>::operator=(disp.value());
        }
    }

    fixedValueFvPatchVectorField::updateCoeffs();

    // If the corresponding point displacement field has a fixedValue type
    // boundary condition, then we wil update it
    setPointDisplacement(*this);
}


Foam::tmp<Foam::Field<Foam::vector>>
Foam::fixedDisplacementFvPatchVectorField::snGrad() const
{
    //- fixedValue snGrad with no correction
    //  return (*this - patchInternalField())*this->patch().deltaCoeffs();

    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + internalField().name() + ")"
        );

    // Unit normals
    vectorField n(patch().nf());

    // Delta vectors
    vectorField delta(patch().delta());

    //- Non-orthogonal correction vectors
    vectorField k(((I - sqr(n)) & delta));

    return
    (
        *this
        - (patchInternalField() + (k & gradField.patchInternalField()))
    )*patch().deltaCoeffs();
}

Foam::tmp<Foam::Field<Foam::vector>>
Foam::fixedDisplacementFvPatchVectorField::gradientBoundaryCoeffs() const
{
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + internalField().name() + ")"
        );

    // Unit normals
    vectorField n(patch().nf());

    // Delta vectors
    vectorField delta(patch().delta());

    //- Non-orthogonal correction vectors
    vectorField k(((I - sqr(n)) & delta));

    return
    (
        patch().deltaCoeffs()
       *(*this - (k & gradField.patchInternalField()))
    );
}

void Foam::fixedDisplacementFvPatchVectorField::write(Ostream& os) const
{
    if (dispSeries_.valid())
    {
        writeEntry(os, dispSeries_());
    }

    fixedValueFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        fixedDisplacementFvPatchVectorField
    );
}


// ************************************************************************* //
