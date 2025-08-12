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
#include "crackerFvMeshTopoChanger.H"
#include "generalFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

template<class SolidTraction>
void SolidCohesiveFvPatchVectorField<SolidTraction>::updateDelta()
{
    // Take a copy of previous delta for under-relaxation
    const vectorField prevDelta = delta_;

    // Take a copy of patch displacement
    vectorField disp = *this;

    const fvMesh& mesh = this->patch().boundaryMesh().mesh();

    // Add on accumulated old displacement if incremental
    if (this->internalField().name() == "DD")
    {
        disp +=
            mesh.lookupObject<volVectorField>
            (
                "D"
            ).oldTime().boundaryField()[this->patch().index()];
    }

    // Cast mesh to a crackerMesh

    if
    (
        !mesh.foundObject<fvMeshTopoChangers::cracker>
        (
            fvMeshTopoChangers::cracker::typeName
        )
    )
    {
        FatalErrorInFunction
            << "Must use " << fvMeshTopoChangers::cracker::typeName
            << " topoChanger" << endl
            << abort(FatalError);
    }

    const fvMeshTopoChangers::cracker& cracker =
        mesh.lookupObject<fvMeshTopoChangers::cracker>
        (
            fvMeshTopoChangers::cracker::typeName
        );

    // Get global crack patch displacement field
    const vectorField globalDisp(cracker.globalCrackField(disp));

    // Update delta
    const labelList& gcfa = cracker.globalCrackFaceAddressing();
    label globalIndex = cracker.localCrackStart();

    forAll(delta_, faceI)
    {
        delta_[faceI] =
            globalDisp[gcfa[globalIndex]]
          - globalDisp[globalIndex];

        globalIndex++;
    }

    // Under-relaxation
    delta_ =
        relaxationFactorDelta_*delta_
      + (1.0 - relaxationFactorDelta_)*prevDelta;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SolidTraction>
SolidCohesiveFvPatchVectorField<SolidTraction>::
SolidCohesiveFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    SolidTraction(p, iF),
    cohesiveZoneModelMaster(p),
    delta_(p.size(), vector::zero),
    relaxationFactorDelta_(1.0),
    relaxationFactorTraction_(1.0),
    curTimeIndex_(-1)
{}


template<class SolidTraction>
SolidCohesiveFvPatchVectorField<SolidTraction>::
SolidCohesiveFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    SolidTraction(p, iF, dict),
    cohesiveZoneModelMaster(p, dict),
    delta_(p.size(), vector::zero),
    relaxationFactorDelta_
    (
        dict.lookupOrDefault<scalar>("relaxationFactorDelta", 1.0)
    ),
    relaxationFactorTraction_
    (
        dict.lookupOrDefault<scalar>("relaxationFactorTraction", 1.0)
    ),
    curTimeIndex_(-1)
{
    if (dict.found("delta"))
    {
        delta_ = vectorField("delta", dict, p.size());
    }
}


template<class SolidTraction>
SolidCohesiveFvPatchVectorField<SolidTraction>::
SolidCohesiveFvPatchVectorField
(
    const SolidCohesiveFvPatchVectorField& cpf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:
    SolidTraction(cpf, p, iF, mapper),
    cohesiveZoneModelMaster(p, cpf),
    delta_(mapper(cpf.delta_)),
    relaxationFactorDelta_(cpf.relaxationFactorDelta_),
    relaxationFactorTraction_(cpf.relaxationFactorTraction_),
    curTimeIndex_(cpf.curTimeIndex_)
{}


template<class SolidTraction>
SolidCohesiveFvPatchVectorField<SolidTraction>::
SolidCohesiveFvPatchVectorField
(
    const SolidCohesiveFvPatchVectorField& cpf,
    const DimensionedField<vector, volMesh>& iF
)
:
    SolidTraction(cpf, iF),
    cohesiveZoneModelMaster(this->patch(), cpf),
    delta_(cpf.delta_),
    relaxationFactorDelta_(cpf.relaxationFactorDelta_),
    relaxationFactorTraction_(cpf.relaxationFactorTraction_),
    curTimeIndex_(cpf.curTimeIndex_)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SolidTraction>
void SolidCohesiveFvPatchVectorField<SolidTraction>::map
(
    const fvPatchVectorField& ptf,
    const fieldMapper& mapper
)
{
    SolidTraction::map(ptf, mapper);

    const SolidCohesiveFvPatchVectorField<SolidTraction>& scfpvf =
        dynamicCast<const SolidCohesiveFvPatchVectorField<SolidTraction>>(ptf);

    const label nOldFaces = delta_.size();
    const label nNewFaces = this->traction().size() - nOldFaces;

    cohesiveZoneModelMaster::map(scfpvf, mapper);
    mapper(delta_, delta_);

    // Only perform mapping if the number of faces on the patch has changed

    vectorField& traction = this->traction();
    scalarField& pressure = this->pressure();
    if
    (
        nNewFaces > 0
     && isA<generalFieldMapper>(mapper)
     && dynamicCast<const generalFieldMapper&>(mapper).direct()
    )
    {

        const labelList& addressing =
            dynamicCast<const generalFieldMapper&>
            (
                mapper
            ).directAddressing();
        const label patchSize = this->patch().size();

        if (patchSize == 1 && nNewFaces == 1)
        {
            label i = 0;

            traction[i] = vector::zero;
            pressure[i] = 0.0;
            delta_[i] = vector::zero;
        }
        else if (patchSize == 2 && nNewFaces == 1)
        {
            label i = 1;

            traction[i] = vector::zero;
            pressure[i] = 0.0;
            delta_[i] = vector::zero;
        }
        else if (patchSize == 2 && nNewFaces == 2)
        {
            label i = 0;

            traction[i] = vector::zero;
            pressure[i] = 0.0;
            delta_[i] = vector::zero;

            i = 1;

            traction[i] = vector::zero;
            pressure[i] = 0.0;
            delta_[i] = vector::zero;
        }
        else
        {
            for (label i = 1; i < patchSize; i++)
            {
                if
                (
                    addressing[i] == 0
                 || (addressing[i] < 0 && i >= nOldFaces)
                )
                {
                    traction[i] = vector::zero;
                    pressure[i] = 0.0;
                    delta_[i] = vector::zero;
                }
            }
        }
    }
}


template<class SolidTraction>
void SolidCohesiveFvPatchVectorField<SolidTraction>::reset
(
    const fvPatchVectorField& ptf
)
{
    SolidTraction::reset(ptf);

    const SolidCohesiveFvPatchVectorField<SolidTraction>& scfpvf =
        refCast<const SolidCohesiveFvPatchVectorField<SolidTraction>>(ptf);

    cohesiveZoneModelMaster::reset(scfpvf);
    delta_.reset(scfpvf.delta_);
}


template<class SolidTraction>
bool SolidCohesiveFvPatchVectorField<SolidTraction>::updateFields()
{
    // Check if it is a new time-step
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        curTimeIndex_ = this->db().time().timeIndex();

        // Update old values within the cohesive zone
        cohesiveZone().updateOldFields();
    }

    const vectorField prevTraction(this->traction());

    bool updateTraction = SolidTraction::updateFields();

    // Update deltas
    updateDelta();

    // Update energies
    cohesiveZone().updateEnergy(this->traction(), delta_);

    // Update and relax tractions
    cohesiveZone().updateTraction
    (
        this->traction(),
        delta_,
        updateTraction
    );

    if (relaxationFactorTraction_ != 1)
    {
        this->traction() =
            relaxationFactorTraction_*this->traction()
            + (1.0 - relaxationFactorTraction_)*prevTraction;
    }

    this->updateForce();

    return true;
}


template<class SolidTraction>
void SolidCohesiveFvPatchVectorField<SolidTraction>::write(Ostream& os) const
{
    SolidTraction::write(os);
    cohesiveZoneModelMaster::write(os);

    writeEntry(os, "delta", delta_);

    writeEntry(os, "relaxationFactorDelta", relaxationFactorDelta_);
    writeEntry(os, "relaxationFactorTraction", relaxationFactorTraction_);
    writeEntry(os, "curTimeIndex", curTimeIndex_);

    writeEntry(os, "value", *this);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
