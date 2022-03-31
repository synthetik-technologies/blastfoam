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

#include "solidTractionRotationFvPatchVectorField.H"
#include "lookupSolidModel.H"
#include "RodriguesRotation.H"
#include "fixedValuePointPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidTractionRotationFvPatchVectorField::
solidTractionRotationFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(p, iF),
    rotationOrigin_(Zero),
    rotationAxis_(Zero),
    origFaceCentres_(p.patch().faceCentres()),
    origPatchPoints_(p.patch().localPoints())
{}


solidTractionRotationFvPatchVectorField::
solidTractionRotationFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidTractionFvPatchVectorField(p, iF),
    rotationOrigin_(dict.lookup<vector>("rotationOrigin")),
    rotationAxis_(dict.lookup<vector>("rotationAxis")),
    origFaceCentres_
    (
        dict.lookupOrDefault<vectorField>
        (
            "origFaceCentres",
            p.patch().faceCentres()
        )
    ),
    origPatchPoints_(p.patch().localPoints())
{
    rotationAxis_ /= mag(rotationAxis_);
    origFaceCentres_ -= (rotationAxis_ & origFaceCentres_)*rotationAxis_;
    origPatchPoints_ -= (rotationAxis_ & origPatchPoints_)*rotationAxis_;
}


solidTractionRotationFvPatchVectorField::
solidTractionRotationFvPatchVectorField
(
    const solidTractionRotationFvPatchVectorField& stpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidTractionFvPatchVectorField(stpvf, p, iF, mapper),
    rotationOrigin_(stpvf.rotationOrigin_),
    rotationAxis_(stpvf.rotationAxis_),
    origFaceCentres_(mapper(stpvf.origFaceCentres_)),
    origPatchPoints_(stpvf.origPatchPoints_)
{}


solidTractionRotationFvPatchVectorField::
solidTractionRotationFvPatchVectorField
(
    const solidTractionRotationFvPatchVectorField& stpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(stpvf, iF),
    rotationOrigin_(stpvf.rotationOrigin_),
    rotationAxis_(stpvf.rotationAxis_),
    origFaceCentres_(stpvf.origFaceCentres_),
    origPatchPoints_(stpvf.origPatchPoints_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidTractionRotationFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    solidTractionFvPatchVectorField::autoMap(m);
    m(origFaceCentres_, origFaceCentres_);
}


void solidTractionRotationFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    solidTractionFvPatchVectorField::rmap(ptf, addr);

    const solidTractionRotationFvPatchVectorField& dmptf =
        refCast<const solidTractionRotationFvPatchVectorField>(ptf);

    origFaceCentres_.rmap(dmptf.origFaceCentres_, addr);
}


void solidTractionRotationFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    solidTractionFvPatchVectorField::updateCoeffs();
    gradient() *=
        mag
        (
            origFaceCentres_ - rotationOrigin_
        );
//     const fvMesh& mesh = patch().boundaryMesh().mesh();
//     const volVectorField& D = mesh.lookupObject<volVectorField>("D");
//     const vectorField& pDOld(D.oldTime().boundaryField()[patch().index()]);
//     vectorField relOrigFaceCentres(origFaceCentres_ - rotationOrigin_);
//
//     bool moving = lookupSolidModel(mesh).movingMesh();
//
//     // Calculate the displacement
//
//     // Lookup the gradient field
//     const fvPatchField<tensor>& gradField =
//         patch().lookupPatchField<volTensorField, tensor>
//         (
//             "grad(" + internalField().name() + ")"
//         );
//     // Lookup the solidModel object
//     const solidModel& solMod =
//         lookupSolidModel(patch().boundaryMesh().mesh());
//
//     // Face unit normals
//     const vectorField n(patch().nf());
//
//     // Delta vectors
//     const vectorField delta(patch().delta());
//
//     // Non-orthogonal correction vectors
//     const vectorField k((I - sqr(n)) & delta);
//
//     vectorField pDD
//     (
//         this->patchInternalField()
//       + (k & gradField.patchInternalField())
//       + solMod.tractionBoundarySnGrad
//         (
//             vectorField(k.size(), Zero),
//             scalarField(k.size(), 0.0),
//             patch()
//         )/patch().deltaCoeffs()
//     );
//
//     tmp<vectorField> newPointsTmp;
//
//     // Calculate total displacement
//     if (internalField().name() == "DD")
//     {
//         pDD += pDOld;
//     }
//
//     if (moving)
//     {
//         newPointsTmp =
//             patch().patch().faceCentres() - rotationOrigin_ + pDD;
//     }
//     else
//     {
//         newPointsTmp = relOrigFaceCentres + pDD;
//     }
//     vectorField& newPoints = newPointsTmp.ref();
//
//     //- Remove displacement in the axis direction
//     newPoints -= (rotationAxis_ & newPoints)*rotationAxis_;
//
//
//     const scalarField& magSf(patch().magSf());
//     const scalar sumMagSf(gSum(magSf));
//     vectorField cross((relOrigFaceCentres ^ newPoints));
//     scalarField c(mag(cross)*(sign(cross & rotationAxis_)));
//     scalarField dot(newPoints & relOrigFaceCentres);
//     scalarField angles(atan2(c, dot));
//
//     scalar sinAngle(gSum(sin(angles)*magSf)/sumMagSf);
//     scalar cosAngle(gSum(cos(angles)*magSf)/sumMagSf);
//     scalar angle(maxMagSqr(angles));//atan2(sinAngle, cosAngle));
//     tensor R(RodriguesRotation(rotationAxis_, angle, false));
//
//     vectorField newFaceCentres
//     (
//         (R & (relOrigFaceCentres)) + rotationOrigin_
//     );
//
//     vectorField disp(newFaceCentres - origFaceCentres_);
//     if (internalField().name() == "DD")
//     {
//         disp -= D.oldTime().boundaryField()[patch().index()];
//     }
//
//     // Set the face displacement
//     fvPatchField<vector>::operator=(disp);
//
//     setPointDisplacement(disp);
//
//     fixedValueFvPatchVectorField::updateCoeffs();
}


void solidTractionRotationFvPatchVectorField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const volVectorField& D = mesh.lookupObject<volVectorField>("D");
    const vectorField& pDOld(D.oldTime().boundaryField()[patch().index()]);
    vectorField relOrigFaceCentres(origFaceCentres_ - rotationOrigin_);

    bool moving = lookupSolidModel(mesh).movingMesh();

    solidTractionFvPatchVectorField::evaluate(commsType);

    // Calculate the displacement
    vectorField pD(*this);
    tmp<vectorField> newPointsTmp;

    // Calculate total displacement
    if (internalField().name() == "DD")
    {
        pD += pDOld;
    }

    if (moving)
    {
        newPointsTmp =
            patch().patch().faceCentres() - rotationOrigin_ + pD;
    }
    else
    {
        newPointsTmp = relOrigFaceCentres + pD;
    }
    vectorField& newPoints = newPointsTmp.ref();

    //- Remove displacement in the axis direction
    newPoints -= (rotationAxis_ & newPoints)*rotationAxis_;


    const scalarField& magSf(patch().magSf());
    const scalar sumMagSf(gSum(magSf));
    vectorField cross((relOrigFaceCentres ^ newPoints));
    scalarField c(mag(cross)*(sign(cross & rotationAxis_)));
    scalarField dot(newPoints & relOrigFaceCentres);
    scalarField angles(atan2(c, dot));

    scalar sinAngle(gSum(sin(angles)*magSf)/sumMagSf);
    scalar cosAngle(gSum(cos(angles)*magSf)/sumMagSf);
    scalar angle(maxMagSqr(angles));//atan2(sinAngle, cosAngle));
    tensor R(RodriguesRotation(rotationAxis_, angle, false));

    vectorField newFaceCentres
    (
        (R & (relOrigFaceCentres)) + rotationOrigin_
    );

    vectorField disp(newFaceCentres - origFaceCentres_);
    if (internalField().name() == "DD")
    {
        disp -= D.oldTime().boundaryField()[patch().index()];
    }

    // Set the face displacement
    fvPatchField<vector>::operator=(disp);

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
                (R & (origPatchPoints_ - rotationOrigin_))
              + rotationOrigin_
            );

            vectorField pointDisp(newPatchPoints - origPatchPoints_);

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
}


void solidTractionRotationFvPatchVectorField::write(Ostream& os) const
{
    solidTractionFvPatchVectorField::write(os);
//     fixedDisplacementFvPatchVectorField::write(os);

    writeEntry(os, "rotationOrigin", rotationOrigin_);
    writeEntry(os, "rotationAxis", rotationAxis_);
    writeEntry(os, "origFaceCentres", origFaceCentres_);
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    solidTractionRotationFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
