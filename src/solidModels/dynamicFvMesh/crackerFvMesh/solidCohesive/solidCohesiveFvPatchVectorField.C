/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "solidCohesiveFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "lookupSolidModel.H"
#include "directFvPatchFieldMapper.H"
#include "crackerFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //


void solidCohesiveFvPatchVectorField::updateDelta()
{
    // Take a copy of previous delta for under-relaxation
    const vectorField prevDelta = delta_;

    // Take a copy of patch displacement
    vectorField disp = *this;

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    // Add on accumulated old displacement if incremental
    if (incremental_)
    {
        const volVectorField& D_0 =
            mesh.objectRegistry::lookupObject<volVectorField>("D").oldTime();

        disp += D_0.boundaryField()[patch().index()];
    }

    // Cast mesh to a crackerMesh

    if (!isA<crackerFvMesh>(mesh))
    {
        FatalErrorInFunction
            << "Mesh should be of type: " << crackerFvMesh::typeName
            << abort(FatalError);
    }

    const crackerFvMesh& crackerMesh =
        dynamicCast<const crackerFvMesh>(mesh);

    // Get global crack patch displacement field
    const vectorField globalDisp(crackerMesh.globalCrackField(disp));

    // Update delta
    const labelList& gcfa = crackerMesh.globalCrackFaceAddressing();
    label globalIndex = crackerMesh.localCrackStart();
    forAll(delta_, faceI)
    {
        delta_[faceI] =
            globalDisp[gcfa[globalIndex]]
          - globalDisp[globalIndex];

        globalIndex++;
    }

    // Under-relaxation
    delta_ = relaxationFactor_*delta_ + (1.0 - relaxationFactor_)*prevDelta;
}


void solidCohesiveFvPatchVectorField::calcCohesiveZone() const
{
    if (cohesiveZoneModelPtr_.valid())
    {
        FatalErrorIn
        (
            "void solidCohesiveFvPatchVectorField::calcCohesiveZone() const"
        )   << "pointer already set" << abort(FatalError);
    }

    if (fieldName_ != "D" && fieldName_ != "DD")
    {
        FatalErrorIn
        (
            "void solidCohesiveFvPatchVectorField::calcCohesiveZone() const"
        )   << "fieldName is not D or DD" << abort(FatalError);
    }

    cohesiveZoneModelPtr_ =
        cohesiveZoneModel::New
        (
            "type", patch(), dict_.subDict("cohesiveZoneModel")
        );
}


const cohesiveZoneModel& solidCohesiveFvPatchVectorField::cohesiveZone() const
{
    if (!cohesiveZoneModelPtr_.valid())
    {
        calcCohesiveZone();
    }

    return cohesiveZoneModelPtr_();
}


cohesiveZoneModel& solidCohesiveFvPatchVectorField::cohesiveZone()
{
    if (!cohesiveZoneModelPtr_.valid())
    {
        calcCohesiveZone();
    }

    return cohesiveZoneModelPtr_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidCohesiveFvPatchVectorField::solidCohesiveFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    fieldName_(internalField().name()),
    cohesiveZoneModelPtr_(NULL),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    delta_(p.size(), vector::zero),
    relaxationFactor_(1.0),
    relaxationFactorTrac_(1.0),
    curTimeIndex_(-1),
    incremental_(false),
    dict_(NULL)
{}


solidCohesiveFvPatchVectorField::solidCohesiveFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    fieldName_(internalField().name()),
    cohesiveZoneModelPtr_(NULL),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    delta_(p.size(), vector::zero),
    relaxationFactor_(dict.lookupOrDefault<scalar>("relaxationFactor", 0.05)),
    relaxationFactorTrac_
    (
        dict.lookupOrDefault<scalar>("relaxationFactorTraction", 0.5)
    ),
    curTimeIndex_(-1),
    incremental_(fieldName_ == "DD"),
    dict_(dict)
{
    if (dict.found("gradient"))
    {
        gradient() = vectorField("gradient", dict, p.size());
    }
    else
    {
        gradient() = vector::zero;
    }

    if (dict.found("value"))
    {
        vectorField::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        vectorField::operator=(vector::zero);
    }

    if (dict.found("traction"))
    {
        traction_ =
            vectorField("traction", dict, p.size());
    }

    if (dict.found("pressure"))
    {
        pressure_ =
            scalarField("pressure", dict, p.size());
    }

    if (dict.found("delta"))
    {
        delta_ = vectorField("delta", dict, p.size());
    }
}


solidCohesiveFvPatchVectorField::solidCohesiveFvPatchVectorField
(
    const solidCohesiveFvPatchVectorField& cpf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(cpf, p, iF, mapper),
    fieldName_(cpf.fieldName_),
    cohesiveZoneModelPtr_(NULL),
    traction_(mapper(cpf.traction_)),
    pressure_(mapper(cpf.pressure_)),
    delta_(mapper(cpf.delta_)),
    relaxationFactor_(cpf.relaxationFactor_),
    relaxationFactorTrac_(cpf.relaxationFactorTrac_),
    curTimeIndex_(cpf.curTimeIndex_),
    incremental_(cpf.incremental_),
    dict_(cpf.dict_)
{
    if (cpf.cohesiveZoneModelPtr_.valid())
    {
        cohesiveZoneModelPtr_.set(cpf.cohesiveZoneModelPtr_->clone().ptr());
    }
}


solidCohesiveFvPatchVectorField::solidCohesiveFvPatchVectorField
(
    const solidCohesiveFvPatchVectorField& cpf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(cpf, iF),
    fieldName_(cpf.fieldName_),
    cohesiveZoneModelPtr_(NULL),
    traction_(cpf.traction_),
    pressure_(cpf.pressure_),
    delta_(cpf.delta_),
    relaxationFactor_(cpf.relaxationFactor_),
    relaxationFactorTrac_(cpf.relaxationFactorTrac_),
    curTimeIndex_(cpf.curTimeIndex_),
    incremental_(cpf.incremental_),
    dict_(cpf.dict_)
{
    if (cpf.cohesiveZoneModelPtr_.valid())
    {
        cohesiveZoneModelPtr_.set(cpf.cohesiveZoneModelPtr_->clone().ptr());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void solidCohesiveFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);

    if (cohesiveZoneModelPtr_.valid())
    {
        cohesiveZone().autoMap(m);
    }

    m(traction_, traction_);
    m(pressure_, pressure_);

    const label nNewFaces = traction_.size() - delta_.size();

    m(delta_, delta_);

    // Only perform mapping if the number of faces on the patch has changed

    if (nNewFaces > 0 && isA<directFvPatchFieldMapper>(m))
    {
        // Reset values on new faces to zero
        // Note: the method below is used to find which faces are new on the
        // patch

        const directFvPatchFieldMapper& dm =
            dynamicCast<const directFvPatchFieldMapper&>(m);
        const labelList& addressing = dm.addressing();
        const label patchSize = patch().size();

        if (patchSize == 1 && nNewFaces == 1)
        {
            label i = 0;

            traction_[i] = vector::zero;
            pressure_[i] = 0.0;
            delta_[i] = vector::zero;
        }
        else if (patchSize == 2 && nNewFaces == 1)
        {
            label i = 1;

            traction_[i] = vector::zero;
            pressure_[i] = 0.0;
            delta_[i] = vector::zero;
        }
        else if (patchSize == 2 && nNewFaces == 2)
        {
            label i = 0;

            traction_[i] = vector::zero;
            pressure_[i] = 0.0;
            delta_[i] = vector::zero;

            i = 1;

            traction_[i] = vector::zero;
            pressure_[i] = 0.0;
            delta_[i] = vector::zero;
        }
        else
        {
            for (label i = 1; i < patchSize; i++)
            {
                if (addressing[i] == 0)
                {
                    traction_[i] = vector::zero;
                    pressure_[i] = 0.0;
                    delta_[i] = vector::zero;
                }
            }
        }
    }
}


// Reverse-map the given fvPatchField onto this fvPatchField
void solidCohesiveFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const solidCohesiveFvPatchVectorField& scptf =
        refCast<const solidCohesiveFvPatchVectorField>(ptf);

    if (cohesiveZoneModelPtr_.valid())
    {
        cohesiveZone().rmap(scptf, addr);
    }

    traction_.rmap(scptf.traction_, addr);
    pressure_.rmap(scptf.pressure_, addr);
    delta_.rmap(scptf.delta_, addr);
}


// Update the coefficients associated with the patch field
void solidCohesiveFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (fieldName_ ==  "D" || fieldName_ ==  "DD")
    {
        // Check if it is a new time-step
        if (curTimeIndex_ != this->db().time().timeIndex())
        {
            curTimeIndex_ = this->db().time().timeIndex();

            // Update old values within the cohesive zone
            cohesiveZone().updateOldFields();
        }

        // Update deltas
        updateDelta();

        // Update energies
        cohesiveZone().updateEnergy(traction_, delta_);

        // Update and relax tractions
        const vectorField prevTraction(traction_);

        cohesiveZone().updateTraction(traction_, delta_);

        traction_ =
            relaxationFactorTrac_*traction_
          + (1.0 - relaxationFactorTrac_)*prevTraction;

        // Lookup the solidModel object
        const solidModel& solMod =
            lookupSolidModel(patch().boundaryMesh().mesh());

        // Set the patch tractions
        gradient() =
            solMod.tractionBoundarySnGrad
            (
                traction_,
                pressure_,
                patch()
            );
    }

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void solidCohesiveFvPatchVectorField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // Calculate boundary value using non-orthogonal correction

    // Lookup explicit gradient
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + fieldName_ + ")"
        );

    // Patch unit normals
    const vectorField n(patch().nf());

    // Delta vectors
    const vectorField delta(patch().delta());

    // Correction vectors
    const vectorField k(delta - n*(n&delta));

    // Set value on the boundary
    Field<vector>::operator=
    (
        this->patchInternalField()
      + (k & gradField.patchInternalField())
      + gradient()/this->patch().deltaCoeffs()
    );

    fvPatchField<vector>::evaluate();
}


void solidCohesiveFvPatchVectorField::write(Ostream& os) const
{
    fixedGradientFvPatchVectorField::write(os);

    if (cohesiveZoneModelPtr_.valid())
    {
        os.writeKeyword("cohesiveZoneModel") << nl;
        os << indent << token::BEGIN_BLOCK << nl;
        cohesiveZone().write(os);
        os << indent << token::END_BLOCK << nl;
    }

    writeEntry(os, "traction", traction_);
    writeEntry(os, "pressure", pressure_);
    writeEntry(os, "delta", delta_);

    writeEntry(os, "relaxationFactor", relaxationFactor_);
    writeEntry(os, "relaxationFactorTraction", relaxationFactorTrac_);
    writeEntry(os, "curTimeIndex", curTimeIndex_);

    writeEntry(os, "value", *this);

    // vectorField posVecs = patch().Cf();
    // forAll(traction_, i)
    // {
    //     Info<< i << " " << traction_[i].x() << " " << delta_[i].x()
    //         << " " << posVecs[i].y() << endl;
    // }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, solidCohesiveFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
