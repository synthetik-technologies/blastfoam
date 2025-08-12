/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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

InClass
    simpleCohesiveZoneFvPatchVectorField

\*---------------------------------------------------------------------------*/

#include "simpleCohesiveZoneFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "pointFields.H"
#include "wallFvPatch.H"
#include "lookupSolidModel.H"
#include "simpleCrackerFvMeshTopoChanger.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void simpleCohesiveZoneFvPatchVectorField::makeSimpleCohesiveZoneLaw() const
{
    if (cohesiveLawPtr_.valid())
    {
        FatalErrorInFunction
            << "Pointer already set!" << abort(FatalError);
    }

    cohesiveLawPtr_.set
    (
        simpleCohesiveZoneLaw::New
        (
            dict_.lookup("simpleCohesiveZoneLaw"),
            dict_
        ).ptr()
    );
}


const simpleCohesiveZoneLaw& simpleCohesiveZoneFvPatchVectorField::law() const
{
    if (cohesiveLawPtr_.empty())
    {
        makeSimpleCohesiveZoneLaw();
    }

    return cohesiveLawPtr_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

simpleCohesiveZoneFvPatchVectorField
::simpleCohesiveZoneFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
    tractionBase(p),
    dict_(),
    totRefValue_(p.size(), vector::zero),
    cohesiveLawPtr_(),
    crackIndicator_(p.size(), 0.0),
    damageIndicator_(p.size(), 0.0),
    relaxationFactor_(1.0),
    separationDistance_(p.size(), 0.0),
    oldSeparationDistance_(p.size(), 0.0),
    unloadingSeparationDistance_(p.size(), 0.0),
    explicitSeparationDistance_(false),
    curTimeIndex_(-1),
    initiationTraction_(p.size(), vector::zero),
    breakOnlyOneFace_(false)
{}


simpleCohesiveZoneFvPatchVectorField
::simpleCohesiveZoneFvPatchVectorField
(
    const simpleCohesiveZoneFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidDirectionMixedFvPatchVectorField(ptf, p, iF, mapper),
    tractionBase(ptf, p, mapper),
    dict_(ptf.dict_),
    totRefValue_(ptf.totRefValue_),
    cohesiveLawPtr_(ptf.cohesiveLawPtr_),
    crackIndicator_(ptf.crackIndicator_),
    damageIndicator_(ptf.damageIndicator_),
    relaxationFactor_(ptf.relaxationFactor_),
    separationDistance_(ptf.separationDistance_),
    oldSeparationDistance_(ptf.oldSeparationDistance_),
    unloadingSeparationDistance_(ptf.unloadingSeparationDistance_),
    explicitSeparationDistance_(ptf.explicitSeparationDistance_),
    curTimeIndex_(ptf.curTimeIndex_),
    initiationTraction_(ptf.initiationTraction_),
    breakOnlyOneFace_(ptf.breakOnlyOneFace_)
{}


simpleCohesiveZoneFvPatchVectorField
::simpleCohesiveZoneFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
    tractionBase(p, dict, false),
    dict_(dict),
    totRefValue_(p.size(), vector::zero),
    cohesiveLawPtr_(),
    crackIndicator_(p.size(), 0.0),
    damageIndicator_(p.size(), 0.0),
    relaxationFactor_(readScalar(dict.lookup("relaxationFactor"))),
    separationDistance_(p.size(), 0.0),
    oldSeparationDistance_(p.size(), 0.0),
    unloadingSeparationDistance_(p.size(), 0.0),
    explicitSeparationDistance_
    (
        dict.lookupOrDefault<Switch>("explicitSeparationDistance", false)
    ),
    curTimeIndex_(-1),
    initiationTraction_(p.size(), vector::zero),
    breakOnlyOneFace_(dict.lookupOrDefault<Switch>("breakOnlyOneFace", false))
{
    // if (!isA<simpleCrackerFvMesh>(patch().boundaryMesh().mesh()))
    // {
    //     FatalErrorIn
    //     (
    //         "simpleCohesiveZoneFvPatchVectorField::"
    //         "simpleCohesiveZoneFvPatchVectorField"
    //     )   << "The " << type() << " boundary condition must be used"
    //         << " with the " << simpleCrackerFvMesh::typeName
    //         << " dynamicFvMesh" << abort(FatalError);
    // }

    if (dict.found("totRefValue"))
    {
        totRefValue_ = vectorField("totRefValue", dict, p.size());
    }
    else
    {
        totRefValue_ = vector::zero;
    }

    if (dict.found("refValue"))
    {
        this->refValue() = vectorField("refValue", dict, p.size());
    }
    else
    {
        this->refValue() = vector::zero;
    }

    if (dict.found("refGradient"))
    {
        this->refGrad() = vectorField("refGradient", dict, p.size());
    }
    else
    {
        this->refGrad() = vector::zero;
    }

    if (dict.found("valueFraction"))
    {
        this->valueFraction() =
            symmTensorField("valueFraction", dict, p.size());
    }
    else
    {
        if (patch().type() == wallFvPatch::typeName)
        {
            this->valueFraction() = I;
        }
        else
        {
            // Symmetry plane
            vectorField n(this->patch().nf());
            this->valueFraction() = sqr(n);
        }
    }

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        directionMixedFvPatchVectorField::evaluate();
    }

    if (dict.found("crackIndicator"))
    {
        crackIndicator_ = scalarField("crackIndicator", dict, p.size());
    }

    if (dict.found("damageIndicator"))
    {
        damageIndicator_ = scalarField("damageIndicator", dict, p.size());
    }

    if (dict.found("separationDistance"))
    {
        separationDistance_ =
            scalarField("separationDistance", dict, p.size());
    }

    if (dict.found("oldSeparationDistance"))
    {
        separationDistance_ =
            scalarField("oldSeparationDistance", dict, p.size());
    }
}


simpleCohesiveZoneFvPatchVectorField
::simpleCohesiveZoneFvPatchVectorField
(
    const simpleCohesiveZoneFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(ptf, iF),
    tractionBase(ptf, this->patch()),
    dict_(ptf.dict_),
    totRefValue_(ptf.totRefValue_),
    cohesiveLawPtr_(ptf.cohesiveLawPtr_),
    crackIndicator_(ptf.crackIndicator_),
    damageIndicator_(ptf.damageIndicator_),
    relaxationFactor_(ptf.relaxationFactor_),
    separationDistance_(ptf.separationDistance_),
    oldSeparationDistance_(ptf.oldSeparationDistance_),
    unloadingSeparationDistance_(ptf.unloadingSeparationDistance_),
    explicitSeparationDistance_(ptf.explicitSeparationDistance_),
    curTimeIndex_(ptf.curTimeIndex_),
    initiationTraction_(ptf.initiationTraction_),
    breakOnlyOneFace_(ptf.breakOnlyOneFace_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void simpleCohesiveZoneFvPatchVectorField::map
(
    const fvPatchVectorField& ptf,
    const fieldMapper& mapper
)
{
     solidDirectionMixedFvPatchVectorField::map(ptf, mapper);

    const simpleCohesiveZoneFvPatchVectorField& sczfvpf =
        refCast<const simpleCohesiveZoneFvPatchVectorField>(ptf);

    mapper(totRefValue_, sczfvpf.totRefValue_);
    mapper(crackIndicator_, sczfvpf.crackIndicator_);
    mapper(damageIndicator_, sczfvpf.damageIndicator_);
    mapper(separationDistance_, sczfvpf.separationDistance_);
    mapper(oldSeparationDistance_, sczfvpf.oldSeparationDistance_);
    mapper(unloadingSeparationDistance_, sczfvpf.unloadingSeparationDistance_);
    mapper(initiationTraction_, sczfvpf.initiationTraction_);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void simpleCohesiveZoneFvPatchVectorField::reset
(
    const fvPatchField<vector>& ptf
)
{
    solidDirectionMixedFvPatchVectorField::reset(ptf);

    const simpleCohesiveZoneFvPatchVectorField& sczfvpf =
        refCast<const simpleCohesiveZoneFvPatchVectorField>(ptf);
    tractionBase::reset(sczfvpf);

    if (cohesiveLawPtr_.empty() && sczfvpf.cohesiveLawPtr_.valid())
    {
        cohesiveLawPtr_.set(sczfvpf.cohesiveLawPtr_().clone().ptr());
    }

    totRefValue_.reset(sczfvpf.totRefValue_);
    crackIndicator_.reset(sczfvpf.crackIndicator_);
    damageIndicator_.reset(sczfvpf.damageIndicator_);
    relaxationFactor_ = sczfvpf.relaxationFactor_;
    separationDistance_.reset(sczfvpf.separationDistance_);
    oldSeparationDistance_.reset(sczfvpf.oldSeparationDistance_);
    unloadingSeparationDistance_.reset
    (
        sczfvpf.unloadingSeparationDistance_
    );
    explicitSeparationDistance_ = sczfvpf.explicitSeparationDistance_;
    curTimeIndex_ = sczfvpf.curTimeIndex_;
    initiationTraction_.reset(sczfvpf.initiationTraction_);
    breakOnlyOneFace_ = sczfvpf.breakOnlyOneFace_;
}


label simpleCohesiveZoneFvPatchVectorField::updateCrack()
{
    const word DDName = internalField().name();

    // Lookup sigma from the solver
    const fvPatchField<symmTensor>& curSigma =
        patch().lookupPatchField<volSymmTensorField, symmTensor>("sigma");

    // Patch normal
    const vectorField n(this->patch().nf());

    // Current normal traction
    const scalarField curNormalTraction(n & (n & curSigma));

    // Current tangential traction
    const vectorField curTangentialTraction((I - sqr(n)) & (n & curSigma));

    // New traction to be set on the patch
    vectorField newTraction(patch().size(), vector::zero);

    // Material strength from cohesive law
    const scalar sigmaMax = law().sigmaMax().value();
    symmTensorField& valueFrac = this->valueFraction();

    // Select potential faces to break
    SLList<label> facesToBreakList;
    SLList<scalar> facesToBreakTractionList;

    forAll(curNormalTraction, faceI)
    {
        if
        (
            (magSqr(valueFrac[faceI]) > (1.0 - SMALL))
         && (curNormalTraction[faceI] >= sigmaMax)
        )
        {
            facesToBreakList.insert(faceI);
            facesToBreakTractionList.insert(curNormalTraction[faceI]);
        }
    }

    labelList facesToBreak(facesToBreakList);
    const List<scalar> facesToBreakTraction(facesToBreakTractionList);

    label nFacesToBreak = facesToBreak.size();

    const label gNFacesToBreak = returnReduce(nFacesToBreak, sumOp<label>());

    if (breakOnlyOneFace_ && gNFacesToBreak)
    {
        if (nFacesToBreak > 1)
        {
            nFacesToBreak = 1;
        }

        // Select internal face with maximum normal traction
        label faceToBreakIndex = -1;
        scalar faceToBreakTraction = 0;
        forAll(facesToBreakTraction, faceI)
        {
            if (facesToBreakTraction[faceI] > faceToBreakTraction)
            {
                faceToBreakTraction = facesToBreakTraction[faceI];
                faceToBreakIndex = facesToBreak[faceI];
            }
        }

        const scalar gMaxTraction =
            returnReduce(faceToBreakTraction, maxOp<scalar>());

        if (Pstream::parRun())
        {
            bool procHasFaceToBreak = false;
            if (nFacesToBreak > 0)
            {
                if (mag(gMaxTraction - faceToBreakTraction) < SMALL)
                {
                    // Maximum traction is on this processor
                    procHasFaceToBreak = true;
                }
                else
                {
                    nFacesToBreak = 0;
                }
            }

            // Check if maximum is present on more then one processors

            label procID = Pstream::nProcs();
            if (procHasFaceToBreak)
            {
                procID = Pstream::myProcNo();
            }

            label minProcID =
                returnReduce<label>(procID, minOp<label>());

            if (procID != minProcID)
            {
                nFacesToBreak = 0;
            }
        }

        facesToBreak.setSize(nFacesToBreak);

        if (nFacesToBreak)
        {
            facesToBreak[0] = faceToBreakIndex;
        }
    }

    forAll(facesToBreak, fI)
    {
        const label faceID = facesToBreak[fI];

        DebugInfo
            << "Switching valueFraction to zero for face "
            << faceID << endl;

        // Switch to full traction boundary condition
        valueFrac[faceID] = symmTensor::zero;
        damageIndicator_[faceID] = 1;
        crackIndicator_[faceID] = 0;

        initiationTraction_[faceID] = n[faceID]*law().sigmaMax().value();

        newTraction[faceID] = initiationTraction_[faceID];

        if (debug)
        {
            Pout<< "Crack has started at face: " << faceID << nl
                << "    normal traction: "
                << curNormalTraction[faceID] << nl
                << "    tangential traction: "
                << curTangentialTraction[faceID] << nl
                << "    initiation traction: "
                << initiationTraction_[faceID]
                << endl;
        }
    }

    this->traction() = newTraction;

    // Lookup the solidModel object
    const solidModel& solMod =
        lookupSolidModel(this->patch().boundaryMesh().mesh());

    // Set traction on the patch
    this->refGrad() =
        solMod.tractionBoundarySnGrad
        (
            this->traction(),
            this->pressure(),
            this->patch()
        );

    return returnReduce(nFacesToBreak, sumOp<label>());
}


bool simpleCohesiveZoneFvPatchVectorField::updateFields()
{
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        const scalar deltaC = law().deltaC().value();

        forAll(unloadingSeparationDistance_, faceI)
        {
            scalar curSepDist = separationDistance_[faceI];

            if (explicitSeparationDistance_)
            {
                curSepDist = oldSeparationDistance_[faceI];
            }

            if (curSepDist < 0)
            {
                curSepDist = 0;
            }

            if (curSepDist > deltaC)
            {
                if
                (
                    curSepDist > unloadingSeparationDistance_[faceI]
                )
                {
                    unloadingSeparationDistance_[faceI] = curSepDist;
                }
            }
        }

        oldSeparationDistance_ = separationDistance_;

        totRefValue_ += this->refValue();

        curTimeIndex_ = this->db().time().timeIndex();
    }

    // Lookup sigma from the solver
    const fvPatchField<symmTensor>& curSigma =
        patch().lookupPatchField<volSymmTensorField, symmTensor>("sigma");


    // Calculate current patch displacement
    vectorField curD(patch().size(), vector::zero);
    if (internalField().name() == "DD")
    {
        const volVectorField& Dold =
            patch().boundaryMesh().mesh().lookupObject<volVectorField>
            (
                "D"
            ).oldTime();

        curD = Dold + (*this);
    }
    else
    {
        curD = (*this);
    }

    // Patch normal
    const vectorField n(this->patch().nf());

    // Current normal traction
    const scalarField curNormalTraction(n & (n & curSigma));

    // New traction to be set
    vectorField newTraction(this->patch().size(), vector::zero);

    // Separation distance
    if (explicitSeparationDistance_)
    {
        separationDistance_ = (n & (totRefValue_ - curD));
    }
    else
    {
        scalarField newSeparationDistance(n & (totRefValue_ - curD));

        // If this is a symmetryPlane type boundary, then we double the
        // displacement; this is because half of the energy is dissipated on the
        // other side of the symmetryPlane
        if (isA<wallFvPatch>(this->patch()))
        {
            // Double the delta for symmetryPlane patches
            newSeparationDistance += newSeparationDistance;
        }

        separationDistance_ =
            separationDistance_
          + relaxationFactor_*(newSeparationDistance - separationDistance_);
    }

    // Check if the separation has become unreasonable
    if
    (
        max(mag(separationDistance_))
      > 1e6*sqrt(sum(this->patch().magSf()))
    )
    {
        FatalErrorInFunction
            << "The separation distance has become very large!"
            << " Possibly a crack has propagated through the entire domain"
            << nl
            << "separationDistance_ is " << separationDistance_ << nl
            << "totRefValue_ is " << totRefValue_ << nl
            << "D is " << curD << nl
            << "n & Dincr" << (n & (totRefValue_ - curD)) << nl
            << abort(FatalError);
    }


    // Check crack propagation

    const symmTensorField& valueFrac = this->valueFraction();
    const scalar deltaC = law().deltaC().value();

    forAll(curNormalTraction, faceI)
    {
        vector cohesiveTraction = vector::zero;

        scalar curSepDist = separationDistance_[faceI];

        if (explicitSeparationDistance_)
        {
            curSepDist = oldSeparationDistance_[faceI];
        }

        if (curSepDist < 0)
        {
            curSepDist = 0;
        }

        if (magSqr(valueFrac[faceI]) < SMALL)
        {
            if
            (
                (curSepDist > deltaC)
             || (curSepDist < (unloadingSeparationDistance_[faceI] - SMALL))
            )
            {
                // Traction free
                cohesiveTraction = vector::zero;
                newTraction[faceI] = cohesiveTraction;

                damageIndicator_[faceI] = 0;
                crackIndicator_[faceI] = 1;
            }
            else
            {
                // Calculate cohesive traction from cohesive zone model
                cohesiveTraction =
                    initiationTraction_[faceI]
                   *law().traction(curSepDist)/law().sigmaMax().value();

                newTraction[faceI] = cohesiveTraction;

                damageIndicator_[faceI] = 1;
                crackIndicator_[faceI] = 0;
            }
        }
    }

    this->traction() = newTraction;
    this->pressure() = Zero;
    return true;
}


void simpleCohesiveZoneFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

   updateFields();

    // Lookup the solidModel object
    const solidModel& solMod =
        lookupSolidModel(this->patch().boundaryMesh().mesh());

    // Set traction on the patch
    refGrad() =
        solMod.tractionBoundarySnGrad
        (
            this->traction(),
            this->pressure(),
            this->patch()
        );

    solidDirectionMixedFvPatchVectorField::updateCoeffs();
}


// Write
void simpleCohesiveZoneFvPatchVectorField::write(Ostream& os) const
{
    solidDirectionMixedFvPatchVectorField::write(os);

    writeEntry(os, "totRefValue", totRefValue_);

    os.writeKeyword("simpleCohesiveZoneLaw") << law().type()
        << token::END_STATEMENT << nl;

    writeEntry(os, "damageIndicator", damageIndicator_);
    writeEntry(os, "crackIndicator", crackIndicator_);
    writeEntry(os, "relaxationFactor", relaxationFactor_);

    law().writeDict(os);

    writeEntry
    (
        os,
        "explicitSeparationDistance",
        explicitSeparationDistance_
    );
    writeEntry(os, "breakOnlyOneFace", breakOnlyOneFace_);
    writeEntry(os, "initiationTraction", initiationTraction_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    simpleCohesiveZoneFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
