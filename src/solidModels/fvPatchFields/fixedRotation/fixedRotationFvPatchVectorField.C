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

#include "fixedRotationFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "RodriguesRotation.H"
#include "fixedValuePointPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedRotationFvPatchVectorField::fixedRotationFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    rotationAngle_(0.0),
    rotationAxis_(vector::zero),
    rotationOrigin_(vector::zero),
    origFaceCentres_(0, vector::zero),
    origPatchPoints_(0, vector::zero),
    angleSeries_(),
    dispSeries_(),
    originSeries_()
{}


fixedRotationFvPatchVectorField::fixedRotationFvPatchVectorField
(
    const fixedRotationFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    rotationAngle_(ptf.rotationAngle_),
    rotationAxis_(ptf.rotationAxis_),
    rotationOrigin_(ptf.rotationOrigin_),
    origFaceCentres_(ptf.origFaceCentres_),
    origPatchPoints_(ptf.origPatchPoints_),
    angleSeries_(ptf.angleSeries_, false),
    dispSeries_(ptf.dispSeries_, false),
    originSeries_(ptf.originSeries_, false)
{}


fixedRotationFvPatchVectorField::fixedRotationFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    rotationAngle_(0.0),
    rotationAxis_(dict.lookup("rotationAxis")),
    rotationOrigin_(vector::zero),
    origFaceCentres_(patch().patch().faceCentres()),
    origPatchPoints_(patch().patch().localPoints()),
    angleSeries_(),
    dispSeries_(),
    originSeries_()
{
    // Check if angle is time-varying
    if (dict.found("rotationAngleSeries"))
    {
        angleSeries_ = Function1<scalar>::New("rotationAngleSeries", dict);
    }
    else
    {
        rotationAngle_ = readScalar(dict.lookup("rotationAngle"));
    }

    // Check if there is a time-varying translation
    if (dict.found("displacementSeries"))
    {
        dispSeries_ = Function1<vector>::New("displacementSeries", dict);
    }

    // Check if origin is time-varying
    if (dict.found("rotationOriginSeries"))
    {
        originSeries_ = Function1<vector>::New("rotationOriginSeries", dict);
    }
    else
    {
        rotationOrigin_ = dict.lookup("rotationOrigin");
    }

    if (dict.found("value"))
    {
        Field<vector>::operator==
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        Field<vector>::operator==
        (
            vectorField(p.size(), vector::zero)
        );
    }

    const word& fieldName = internalField().name();

    if
    (
        fieldName != "D" && fieldName != "DD"
     && fieldName != "D_0" && fieldName != "DD_0"
     && fieldName != "D_0_0" && fieldName != "DD_0_0"
     && fieldName != "D_0_0_0" && fieldName != "DD_0_0_0"
    )
    {
        FatalErrorIn
        (
            "fixedRotationFvPatchVectorField::"
            "fixedRotationFvPatchVectorField(...)"
        )   << "The displacement field should be D or DD"
            << ": current name is " << fieldName
            << abort(FatalError);
    }
}


fixedRotationFvPatchVectorField::fixedRotationFvPatchVectorField
(
    const fixedRotationFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    rotationAngle_(pivpvf.rotationAngle_),
    rotationAxis_(pivpvf.rotationAxis_),
    rotationOrigin_(pivpvf.rotationOrigin_),
    origFaceCentres_(pivpvf.origFaceCentres_),
    origPatchPoints_(pivpvf.origPatchPoints_),
    angleSeries_(pivpvf.angleSeries_, false),
    dispSeries_(pivpvf.dispSeries_, false),
    originSeries_(pivpvf.originSeries_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<vector> > fixedRotationFvPatchVectorField::
snGrad() const
{
    // fixedValue snGrad with no correction
    // return (*this - patchInternalField())*this->patch().deltaCoeffs();

    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + internalField().name() + ")"
        );

    vectorField n(this->patch().nf());
    vectorField delta(this->patch().delta());

    // Correction vectors
    vectorField k(delta - n*(n & delta));

    return
    (
        *this
      - (patchInternalField() + (k & gradField.patchInternalField()))
    )*this->patch().deltaCoeffs();
}


void fixedRotationFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (angleSeries_.valid())
    {
        rotationAngle_ = angleSeries_->value(this->db().time().timeOutputValue());
    }

    if (originSeries_.valid())
    {
        rotationOrigin_ = originSeries_->value(this->db().time().timeOutputValue());
    }

    // Rotation tensor
    const tensor rotMat = RodriguesRotation(rotationAxis_, rotationAngle_);

    vectorField newFaceCentres
    (
        (rotMat & (origFaceCentres_ - rotationOrigin_)) + rotationOrigin_
    );

    vectorField disp(newFaceCentres - origFaceCentres_);

    // Superimposed translation
    if (dispSeries_.valid())
    {
        disp += dispSeries_->value(this->db().time().timeOutputValue());
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if (internalField().name() == "DD")
    {
        const volVectorField& D = mesh.lookupObject<volVectorField>("D");

        disp -= D.oldTime().boundaryField()[patch().index()];
    }

    // Set the face displacement
    fvPatchField<vector>::operator==
    (
        disp
    );

    // If the point displacement field is found and is fixedValue then we will
    // update it using const_cast
    if
    (
        mesh.foundObject<pointVectorField>
        (
            "point" + internalField().name()
        )
    )
    {
        pointVectorField& pointDDField =
            const_cast<pointVectorField&>
            (
                mesh.lookupObject<pointVectorField>
                (
                    "point" + internalField().name()
                )
            );

        if
        (
            pointDDField.boundaryField()[patch().index()].type()
         == fixedValuePointPatchVectorField::typeName
        )
        {
            fixedValuePointPatchVectorField& pointDD =
                refCast<fixedValuePointPatchVectorField>
                (
                    pointDDField.boundaryFieldRef()[patch().index()]
                );

            vectorField newPatchPoints
            (
                (rotMat & (origPatchPoints_ - rotationOrigin_))
              + rotationOrigin_
            );

            vectorField pointDisp(newPatchPoints - origPatchPoints_);

            // Superimposed translation
            if (dispSeries_.valid())
            {
                pointDisp +=
                    dispSeries_->value(this->db().time().timeOutputValue());
            }

            const labelList& meshPoints =
                mesh.boundaryMesh()[patch().index()].meshPoints();

            if (internalField().name() == "DD")
            {
                // Lookup the accumulated total displacement
                const pointVectorField& pointDField =
                    mesh.lookupObject<pointVectorField>("pointD").oldTime();

                forAll(meshPoints, pointI)
                {
                    pointDisp[pointI] -= pointDField[meshPoints[pointI]];
                }
            }

            pointDD == pointDisp;
        }
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


tmp<Field<vector> > fixedRotationFvPatchVectorField::
gradientBoundaryCoeffs() const
{
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + internalField().name() + ")"
        );

    vectorField n(this->patch().nf());
    vectorField delta(this->patch().delta());

    //- correction vector
    vectorField k(delta - n*(n & delta));

    return
        this->patch().deltaCoeffs()
       *(
           *this - (k & gradField.patchInternalField())
       );
}


void fixedRotationFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchVectorField::write(os);

    if (angleSeries_.valid())
    {
        writeEntry(os, "rotationAngleSeries", angleSeries_());
    }
    else
    {
        writeEntry(os, "rotationAngle", rotationAngle_);
    }
    if (dispSeries_.valid())
    {
        writeEntry(os, "displacementSeries", dispSeries_());
    }
    if (originSeries_.valid())
    {
        writeEntry(os, "rotationOriginSeries", originSeries_());
    }
    else
    {
        writeEntry(os, "rotationOrigin", rotationOrigin_);
    }
    writeEntry(os, "rotationAxis", rotationAxis_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    fixedRotationFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
