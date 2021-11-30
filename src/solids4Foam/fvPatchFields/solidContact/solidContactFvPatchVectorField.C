/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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
    solidContactFvPatchVectorField

\*---------------------------------------------------------------------------*/

#include "solidContactFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "Switch.H"
#include "pointFields.H"
#include "polyPatchID.H"
#include "ZoneIDs.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solidContactFvPatchVectorField::movingMesh() const
{
    // Check if the solid model moves the mesh
    return lookupSolidModel(patch().boundaryMesh().mesh()).movingMesh();
}


void Foam::solidContactFvPatchVectorField::makeShadowPatchNames
(
    const dictionary& dict
) const
{
    if (master_)
    {
        // Check if only one shadow patch is specified
        if (dict.found("shadowPatch"))
        {
            Info<< "Reading individual shadowPatch" << endl;

            // Just one shadow patch
            shadowPatchNames_.setSize(1);
            shadowPatchNames_[0] = word(dict.lookup("shadowPatch"));
        }
        else if (dict.found("shadowPatches"))
        {
            Info<< "Reading list of shadowPatches" << endl;

            // Shadow patches defined as a list
            shadowPatchNames_ = wordList(dict.lookup("shadowPatches"));
        }
        else
        {
            FatalErrorInFunction
                << "'shadowPatch' OR 'shadowPatches' should be defined"
                << abort(FatalError);
        }

        // It is an error to defined both shadowPatch and shadowPatches
        if (dict.found("shadowPatch") && dict.found("shadowPatches"))
        {
            FatalErrorInFunction
                << "'shadowPatch' OR 'shadowPatches' should be defined: "
                << "not both!" << abort(FatalError);
        }
    }
    else
    {
        // If this is not the master then we will assume there is only one
        // shadow i.e. the master is the shadow
        shadowPatchNames_.setSize(1);
        shadowPatchNames_[0] = word(dict.lookup("shadowPatch"));
    }
}


void Foam::solidContactFvPatchVectorField::calcShadowPatchIndices() const
{
    if (shadowPatchIndicesPtr_)
    {
        FatalErrorInFunction
            << "pointer already set" << abort(FatalError);
    }

    shadowPatchIndicesPtr_ = new labelList(shadowPatchNames().size(), -1);
    labelList& shadowPatchIndices = *shadowPatchIndicesPtr_;

    forAll(shadowPatchIndices, shadPatchI)
    {
        shadowPatchIndices[shadPatchI] =
            patch().patch().boundaryMesh().findPatchID
            (
                shadowPatchNames_[shadPatchI]
            );

        if (shadowPatchIndices[shadPatchI] == -1)
        {
            FatalErrorInFunction
                << "shadowPatch " << shadowPatchNames_[shadPatchI]
                << " not found!" << abort(FatalError);
        }
    }
}


void Foam::solidContactFvPatchVectorField::makeNormalModels
(
    const dictionary& dict
) const
{
    normalModels_.setSize(shadowPatchNames().size());

    forAll (normalModels_, shadPatchI)
    {
        // Check if only one shadow patch is defined
        const dictionary* contactDictPtr = NULL;
        if (normalModels_.size() == 1)
        {
            if (dict.found("normalContactModel") )
            {
                contactDictPtr = &dict;
            }
            else
            {
                contactDictPtr =
                    &dict.subDict
                    (
                        patch().name() + "_to_"
                      + shadowPatchNames()[shadPatchI] + "_dict"
                    );
            }
        }
        else
        {
            contactDictPtr =
                &dict.subDict
                (
                    patch().name() + "_to_"
                  + shadowPatchNames()[shadPatchI] + "_dict"
                );
        }
        const dictionary& contactDict = *contactDictPtr;

        // Create contact model
        normalModels_.set
        (
            shadPatchI,
            normalContactModel::New
            (
                word(contactDict.lookup("normalContactModel")),
                patch(),
                contactDict,
                patch().index(),                  // master
                shadowPatchIndices()[shadPatchI], // slave
                zone().globalPatch(),
                shadowZones()[shadPatchI].globalPatch()
            ).ptr()
        );
    }

    // Initialise penalty scales to -1
    Info<< "    Initialising stored previous normalPenaltyFactors" << endl;
    normalPenaltyFactors_.setSize(normalModels_.size(), -1);
}


void Foam::solidContactFvPatchVectorField::makeFrictionModels
(
    const dictionary& dict
) const
{
    frictionModels_.setSize(shadowPatchNames().size());

    forAll (frictionModels_, shadPatchI)
    {
        // Check if only one shadow patch is defined
        const dictionary* contactDictPtr = NULL;
        if (frictionModels_.size() == 1)
        {
            if (dict.found("frictionContactModel") )
            {
                contactDictPtr = &dict;
            }
            else
            {
                contactDictPtr =
                    &dict.subDict
                    (
                        patch().name() + "_to_"
                      + shadowPatchNames()[shadPatchI] + "_dict"
                    );
            }
        }
        else
        {
            contactDictPtr =
                &dict.subDict
                (
                    patch().name() + "_to_"
                  + shadowPatchNames()[shadPatchI] + "_dict"
                );
        }
        const dictionary& contactDict = *contactDictPtr;

        // Create contact model
        frictionModels_.set
        (
            shadPatchI,
            frictionContactModel::New
            (
                word(contactDict.lookup("frictionContactModel")),
                patch(),
                contactDict,
                patch().index(),                 // master
                shadowPatchIndices()[shadPatchI] // slave
            ).ptr()
        );
    }
}


void Foam::solidContactFvPatchVectorField::clearOut()
{
    if (debug)
    {
        InfoInFunction
            << patch().name() << " : clearOut" << endl;
    }

    deleteDemandDrivenData(shadowPatchIndicesPtr_);
    deleteDemandDrivenData(zonePtr_);
    shadowZones_.clear();
    zoneToZones_.clear();
    scaleTractionFieldPtr_.clear();
}

Foam::scalarField
Foam::solidContactFvPatchVectorField::scaleTractionField() const
{
    if (scaleTractionFieldPtr_.empty())
    {
        makeScaleTractionField();
    }

    return scaleTractionFieldPtr_();
}


void Foam::solidContactFvPatchVectorField::makeScaleTractionField() const
{
    if (scaleTractionFieldPtr_.valid())
    {
        FatalErrorInFunction
            << "Pointer already set!" << abort(FatalError);
    }

    scaleTractionFieldPtr_.set(new scalarField(patch().size(), 1.0));
    scalarField& scaleTractionField = scaleTractionFieldPtr_();

    // Find all faces on the patch that are adjacent to faces on the
    // downstream patch

    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const scalar scaleFactor =
        readScalar(dict_.lookup("downstreamScaleFactor"));

    // Downstream patch name
    const word patchName = dict_.lookup("downstreamPatchName");
    const label patchID = mesh.boundaryMesh().findPatchID(patchName);
    if (patchID == -1)
    {
        FatalErrorIn(type() + "::makeScaleTractionField()")
            << "Cannot find patch " << patchName << abort(FatalError);
    }

    const unallocLabelList& faceCells = patch().faceCells();
    const cellList& cells = mesh.cells();

    forAll(scaleTractionField, fI)
    {
        const label cellID = faceCells[fI];
        const cell& curCell = cells[cellID];

        forAll(curCell, cfI)
        {
            const label cellFaceID = curCell[cfI];

            if (!mesh.isInternalFace(cellFaceID))
            {
                if (mesh.boundaryMesh().whichPatch(cellFaceID) == patchID)
                {
                    scaleTractionField[fI] = scaleFactor;
                    break;
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidContactFvPatchVectorField::solidContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(p, iF),
    dict_(),
    master_(false),
    writeZoneVTK_(false),
    writePointDistanceFields_(false),
    shadowPatchNames_(),
    shadowPatchIndicesPtr_(NULL),
    rigidMaster_(false),
    normalModels_(),
    frictionModels_(),
    normalPenaltyFactors_(),
    zonePtr_(NULL),
    shadowZones_(),
    zoneToZones_(),
    quickReject_(Foam::ggiInterpolation::AABB),
    regionOfInterestTopCorner_(vector::max),
    regionOfInterestBottomCorner_(vector::min),
    regionOfInterest_(vector::min, vector::max),
    contact_(0),
    contactPerShadow_(),
    scaleFaceTractionsNearDownstreamPatch_(false),
    scaleTractionFieldPtr_(),
    curTimeIndex_(-1)
{}


Foam::solidContactFvPatchVectorField::solidContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:   solidTractionFvPatchVectorField(p, iF),
    dict_(dict),
    master_(dict.lookup("master")),
    writeZoneVTK_(dict.lookupOrDefault<Switch>("writeZoneVTK", false)),
    writePointDistanceFields_
    (
        dict.lookupOrDefault<Switch>("writePointDistanceFields", false)
    ),
    shadowPatchNames_(),
    shadowPatchIndicesPtr_(NULL),
    rigidMaster_(false),
    normalModels_(),
    frictionModels_(),
    normalPenaltyFactors_(),
    zonePtr_(NULL),
    shadowZones_(),
    zoneToZones_(),
    quickReject_
    (
        ggiInterpolation::quickRejectNames_
        [
            dict.lookupOrDefault<word>("quickReject", "AABB")
        ]
    ),
    regionOfInterestTopCorner_
    (
        dict.lookupOrDefault<vector>
        (
            "regionOfInterestTopCorner",
            vector::max
        )
    ),
    regionOfInterestBottomCorner_
    (
        dict.lookupOrDefault<vector>
        (
            "regionOfInterestBottomCorner",
            vector::min
        )
    ),
    regionOfInterest_
    (
        boundBox
        (
            regionOfInterestBottomCorner_,
            regionOfInterestTopCorner_
        )
//        dict.lookupOrDefault<boundBox>
//        (
//            "regionOfInterest",
//            boundBox(vector::min, vector::max)
//        )
    ),
    contact_(patch().size(), 0.0),
    contactPerShadow_(),
    scaleFaceTractionsNearDownstreamPatch_
    (
        dict.lookupOrDefault<Switch>
        (
            "scaleFaceTractionsNearDownstreamPatch",
            Switch(false)
        )
    ),
    scaleTractionFieldPtr_(),
    curTimeIndex_(-1)
{
    if (debug)
    {
        Info<< "Creating " << solidContactFvPatchVectorField::typeName
            << " patch" << endl;
    }

    // Master creates contact laws
    if (master_)
    {
        rigidMaster_ = Switch(dict.lookup("rigidMaster"));

        if (debug)
        {
            Info<< "    writePointDistanceFields: " << writePointDistanceFields_
                << endl;
        }

        if (scaleFaceTractionsNearDownstreamPatch_)
        {
            WarningInFunction
                << "scaleFaceTractionsNearDownstreamPatch can only be applied on"
                << "the slave patch: this option will be ignored for the master "
                << "patch!" << endl;
        }
    }

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
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        Field<vector>::operator=
        (
            patchInternalField() + gradient()/patch().deltaCoeffs()
        );
    }
}


Foam::solidContactFvPatchVectorField::solidContactFvPatchVectorField
(
    const solidContactFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidTractionFvPatchVectorField(ptf, p, iF, mapper),
    dict_(ptf.dict_),
    master_(ptf.master_),
    writeZoneVTK_(ptf.writeZoneVTK_),
    writePointDistanceFields_(ptf.writePointDistanceFields_),
    shadowPatchNames_(ptf.shadowPatchNames_),
    shadowPatchIndicesPtr_(NULL),
    rigidMaster_(ptf.rigidMaster_),
    normalModels_(ptf.normalModels_),
    frictionModels_(ptf.frictionModels_),
    normalPenaltyFactors_(ptf.normalPenaltyFactors_.size(), -1),
    zonePtr_(NULL),
    shadowZones_(),
    zoneToZones_(),
    quickReject_(ptf.quickReject_),
    regionOfInterestTopCorner_(ptf.regionOfInterestTopCorner_),
    regionOfInterestBottomCorner_(ptf.regionOfInterestBottomCorner_),
    regionOfInterest_(ptf.regionOfInterest_),
    contact_(ptf.contact_),
    contactPerShadow_(),
    scaleFaceTractionsNearDownstreamPatch_
    (
        ptf.scaleFaceTractionsNearDownstreamPatch_
    ),
    scaleTractionFieldPtr_(),
    curTimeIndex_(ptf.curTimeIndex_)
{
    // Do not copy pointer objects: they will be re-created.
}


Foam::solidContactFvPatchVectorField::solidContactFvPatchVectorField
(
    const solidContactFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(ptf, iF),
    dict_(ptf.dict_),
    master_(ptf.master_),
    writeZoneVTK_(ptf.writeZoneVTK_),
    writePointDistanceFields_(ptf.writePointDistanceFields_),
    shadowPatchNames_(ptf.shadowPatchNames_),
    shadowPatchIndicesPtr_(NULL),
    rigidMaster_(ptf.rigidMaster_),
    normalModels_(ptf.normalModels_),
    frictionModels_(ptf.frictionModels_),
    normalPenaltyFactors_(ptf.normalPenaltyFactors_.size(), -1),
    zonePtr_(NULL),
    shadowZones_(),
    zoneToZones_(),
    quickReject_(ptf.quickReject_),
    regionOfInterestTopCorner_(ptf.regionOfInterestTopCorner_),
    regionOfInterestBottomCorner_(ptf.regionOfInterestBottomCorner_),
    regionOfInterest_(ptf.regionOfInterest_),
    contact_(ptf.contact_),
    contactPerShadow_(),
    scaleFaceTractionsNearDownstreamPatch_
    (
        ptf.scaleFaceTractionsNearDownstreamPatch_
    ),
    scaleTractionFieldPtr_(),
    curTimeIndex_(ptf.curTimeIndex_)
{
    // Do not copy pointer objects
}


// * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * * //


Foam::solidContactFvPatchVectorField::
~solidContactFvPatchVectorField()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solidContactFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    if (debug)
    {
        InfoInFunction
            << nl << "autoMap: field = " << internalField().name()
            << ", patch = " << patch().name() << endl;
    }

    solidTractionFvPatchVectorField::autoMap(m);

    m(contact_, contact_);
    scaleTractionFieldPtr_.clear();

    if (contactPerShadow_.size())
    {
        forAll(contactPerShadow_, shadI)
        {
            m(contactPerShadow_[shadI], contactPerShadow_[shadI]);
        }
    }

    if (shadowPatchNames_.size() > 0)
    {
        // Let the contact models know about the mapping
        // Be careful, we must pass slave
        // FIX PC 21-Sep-17: move this check inside if (shadowPatchNames ... )
        if (!master_)
        {
            normalModelForThisSlave().autoMap(m);
            frictionModelForThisSlave().autoMap(m);
        }
    }

    // Force all data to be re-created when needed
    clearOut();

    // Reset normal pelanty factors to reinitialise normal models
    normalPenaltyFactors_ =  -1;
}


void Foam::solidContactFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    solidTractionFvPatchVectorField::rmap(ptf, addr);

    const solidContactFvPatchVectorField& dmptf =
        refCast<const solidContactFvPatchVectorField>(ptf);

    // PC, I'm not sure if this "if" check should be here...
    if (shadowPatchNames_.size() > 0)
    {
        contact_.rmap(dmptf.contact_, addr);

        if (dmptf.contactPerShadow_.size())
        {
            // Force contactPerShadow to be initialised
            contactPerShadow();

            forAll(contactPerShadow_, shadI)
            {
                contactPerShadow_[shadI].rmap
                (
                    dmptf.contactPerShadow_[shadI], addr
                );
            }
        }
    }

    scaleTractionFieldPtr_.clear();
}


const Foam::wordList&
Foam::solidContactFvPatchVectorField::shadowPatchNames() const
{
    if (shadowPatchNames_.size() == 0)
    {
        makeShadowPatchNames(dict_);
    }

    return shadowPatchNames_;
}


const Foam::labelList&
Foam::solidContactFvPatchVectorField::shadowPatchIndices() const
{
    if (!shadowPatchIndicesPtr_)
    {
        calcShadowPatchIndices();
    }

    return *shadowPatchIndicesPtr_;
}


const Foam::solidContactFvPatchVectorField&
Foam::solidContactFvPatchVectorField::shadowPatchField() const
{
    if (shadowPatchIndices().size() != 1)
    {
        FatalErrorInFunction
            << "This function can only be called for a patch with 1 shadow "
            << "patch; this patch has " << shadowPatchIndices().size()
            << " shadow patches!" << abort(FatalError);
    }

    const volVectorField& field =
        db().lookupObject<volVectorField>(internalField().name());

    return
        refCast<const solidContactFvPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndices()[0]]
        );
}


Foam::PtrList<Foam::normalContactModel>&
Foam::solidContactFvPatchVectorField::normalModels()
{
    if (master_)
    {
        if (normalModels_.size() == 0)
        {
            makeNormalModels(dict_);
        }

        return normalModels_;
    }
    else
    {
        const volVectorField& field =
            db().lookupObject<volVectorField>
            (
                this->internalField().name()
            );

        solidContactFvPatchVectorField& shadowPatchField =
            const_cast<solidContactFvPatchVectorField&>
            (
                refCast<const solidContactFvPatchVectorField>
                (
                    field.boundaryField()[shadowPatchIndices()[0]]
                )
            );

        return shadowPatchField.normalModels();
    }
}


const Foam::PtrList<Foam::normalContactModel>&
Foam::solidContactFvPatchVectorField::normalModels() const
{
    if (master_)
    {
        if (normalModels_.size() == 0)
        {
            makeNormalModels(dict_);
        }

        return normalModels_;
    }
    else
    {
        const volVectorField& field =
            db().lookupObject<volVectorField>
            (
                this->internalField().name()
            );

        const solidContactFvPatchVectorField& shadowPatchField =
            refCast<const solidContactFvPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndices()[0]]
            );

        return shadowPatchField.normalModels();
    }
}


Foam::PtrList<Foam::frictionContactModel>&
Foam::solidContactFvPatchVectorField::frictionModels()
{
    if (master_)
    {
        if (frictionModels_.size() == 0)
        {
            makeFrictionModels(dict_);
        }

        return frictionModels_;
    }
    else
    {
        const volVectorField& field =
            db().lookupObject<volVectorField>
            (
                this->internalField().name()
            );

        solidContactFvPatchVectorField& shadowPatchField =
            const_cast<solidContactFvPatchVectorField&>
            (
                refCast<const solidContactFvPatchVectorField>
                (
                    field.boundaryField()[shadowPatchIndices()[0]]
                )
            );

        return shadowPatchField.frictionModels();
    }
}


const Foam::PtrList<Foam::frictionContactModel>&
Foam::solidContactFvPatchVectorField::frictionModels() const
{
    if (master_)
    {
        if (frictionModels_.size() == 0)
        {
            makeFrictionModels(dict_);
        }

        return frictionModels_;
    }
    else
    {
        const volVectorField& field =
            db().lookupObject<volVectorField>
            (
                this->internalField().name()
            );

        solidContactFvPatchVectorField& shadowPatchField =
            const_cast<solidContactFvPatchVectorField&>
            (
                refCast<const solidContactFvPatchVectorField>
                (
                    field.boundaryField()[shadowPatchIndices()[0]]
                )
            );

        return shadowPatchField.frictionModels();
    }
}


Foam::normalContactModel&
Foam::solidContactFvPatchVectorField::normalModelForThisSlave()
{
    if (master_)
    {
        FatalErrorInFunction
            << "The master is not allowed to called this fucntion!"
            << abort(FatalError);
    }

    // Lookup the master patch corresponding to the current slave patch
    const volVectorField& field =
        db().lookupObject<volVectorField>
        (
            this->internalField().name()
        );

    if (returnReduce(shadowPatchIndices().size() == 0, maxOp<bool>()))
    {
        FatalErrorInFunction
            << "shadowPatchIndices().size() == 0" << exit(FatalError);
    }

    const solidContactFvPatchVectorField& masterPatchField =
        refCast<const solidContactFvPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndices()[0]]
        );

    // The master may have multiple slaves so we need to find which model
    // corresponds to the current slave patch
    const wordList& shadPatchNames = masterPatchField.shadowPatchNames();
    label masterShadowID = -1;
    forAll(shadPatchNames, shadPatchI)
    {
        if (shadPatchNames[shadPatchI] == patch().name())
        {
            masterShadowID = shadPatchI;
            break;
        }
    }

    if (masterShadowID == -1)
    {
        FatalErrorInFunction
            << "Something went wrong when looking for the shadowPatch"
            << abort(FatalError);
    }

    return normalModels()[masterShadowID];
}


Foam::frictionContactModel&
Foam::solidContactFvPatchVectorField::frictionModelForThisSlave()
{
    if (master_)
    {
        FatalErrorInFunction
            << "The master is not allowed to called this function!"
            << abort(FatalError);
    }

    // Lookup the master patch corresponding to the current slave patch
    const volVectorField& field =
        db().lookupObject<volVectorField>
        (
            this->internalField().name()
        );

    const solidContactFvPatchVectorField& masterPatchField =
        refCast<const solidContactFvPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndices()[0]]
        );

    // The master may have multiple slaves so we need to find which model
    // corresponds to the current slave patch
    const wordList& shadPatchNames = masterPatchField.shadowPatchNames();
    label masterShadowID = -1;
    forAll(shadPatchNames, shadPatchI)
    {
        if (shadPatchNames[shadPatchI] == patch().name())
        {
            masterShadowID = shadPatchI;
            break;
        }
    }

    if (masterShadowID == -1)
    {
        FatalErrorInFunction
            << "Something went wrong when looking for the shadowPatch"
            << abort(FatalError);
    }

    return frictionModels()[masterShadowID];
}


void Foam::solidContactFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        // Update old quantities at the start of a new time-step
        curTimeIndex_ = this->db().time().timeIndex();

        if (master_)
        {
            // Let the contact models know that it is a new time-step, in case
            // they need to update anything
            forAll(shadowPatchNames(), shadPatchI)
            {
                normalModels()[shadPatchI].newTimeStep();
                frictionModels()[shadPatchI].newTimeStep();

                // Force N^2 contact search at least once per time-step
                zoneToZones()[shadPatchI].clearPrevCandidateMasterNeighbors();
            }
        }
    }

    // Move the master and slave zone to the deformed configuration
    moveZonesToDeformedConfiguration();

    // Delete the zone-to-zone interpolator weights as the zones have moved
    const wordList& shadPatchNames = shadowPatchNames();
    forAll(shadPatchNames, shadPatchI)
    {
        zoneToZones()[shadPatchI].movePoints
        (
            tensorField(0), tensorField(0), vectorField(0)
        );
    }

    // Calculate and apply contact forces
    if (master_)
    {
        // Reset the traction to zero as we will accumulate it over all the
        // shadow patches
        traction() = vector::zero;

        forAll(shadPatchNames, shadPatchI)
        {
            // Calculate the slave patch face unit normals as they are used by
            // both the normal and friction models
            const vectorField shadowPatchFaceNormals
            (
                shadowZones()[shadPatchI].globalFaceToPatch
                (
                    shadowZones()[shadPatchI].globalPatch().faceNormals()
                )
            );

            // Interpolate the master displacement increment to the slave patch
            // as it is required by specific normal and friction contact models

            vectorField patchDD(patch().size(), vector::zero);
            vectorField shadowPatchDD
            (
                patch().boundaryMesh()[shadowPatchIndices()[shadPatchI]].size(),
                vector::zero
            );

            if (movingMesh())
            {
                // Updated Lagrangian, we will directly lookup the displacement
                // increment

                const volVectorField& DD =
                    db().lookupObject<volVectorField>("DD");

                patchDD = DD.boundaryField()[patch().index()];
                shadowPatchDD =
                    DD.boundaryField()[shadowPatchIndices()[shadPatchI]];
            }
            else
            {
                // We will lookup the total displacement and old total
                // displacement

                const volVectorField& D =
                    db().lookupObject<volVectorField>("D");

                patchDD =
                    D.boundaryField()[patch().index()]
                  - D.oldTime().boundaryField()[patch().index()];
                shadowPatchDD =
                    D.boundaryField()[shadowPatchIndices()[shadPatchI]]
                  - D.oldTime().boundaryField()
                    [
                        shadowPatchIndices()[shadPatchI]
                    ];
            }

            // Master zone DD
            const vectorField zoneDD(zone().patchFaceToGlobal(patchDD));

            // Master patch DD interpolated to the slave patch
            const vectorField patchDDInterpToShadowPatch
            (
                shadowZones()[shadPatchI].globalFaceToPatch
                (
                    zoneToZones()[shadPatchI].masterToSlave(zoneDD)()
                )
            );

            // Calculate normal contact forces
            // shadowPatchDD is the DU on the shadow patch, whereas
            // patchDDInterpToShadowPatch is the master patch DU interpolated to
            // the shadow; and the difference between these two is the slip (and
            // also the normal component of DU)
            normalModels()[shadPatchI].correct
            (
                shadowPatchFaceNormals,
                shadowZones()[shadPatchI].globalPointToPatch
                (
                    zoneToZones()[shadPatchI].slavePointDistanceToIntersection()
                ),
                // zoneToZones()[shadPatchI],
                shadowPatchDD,
                patchDDInterpToShadowPatch
            );

            // Calculate friction contact forces
            frictionModels()[shadPatchI].correct
            (
                normalModels()[shadPatchI].slavePressure(),
                shadowPatchFaceNormals,
                normalModels()[shadPatchI].areaInContact(),
                shadowPatchDD,
                patchDDInterpToShadowPatch
            );

            if (rigidMaster_)
            {
                // Set to master to traction free to mimic a rigid contact
                traction() = vector::zero;

                // Set contact indicator field
                contactPerShadow()[shadPatchI] = 0.0;
            }
            else
            {
                // Interpolate slave traction to the master
                const vectorField slavePatchTraction
                (
                   - frictionModels()[shadPatchI].slaveTractionForMaster()
                   - normalModels()[shadPatchI].slavePressure()
                );

                const vectorField slaveZoneTraction
                (
                    shadowZones()[shadPatchI].patchFaceToGlobal
                    (
                        slavePatchTraction
                    )
                );

                // We have two options for interpolating from the slave to the
                // master:
                // 1. face-to-face
                // 2. point-to-point
                // We will use 1.

                // Calculate traction for this contact
                vectorField tractionForThisShadow
                (
                    zone().globalFaceToPatch
                    (
                        zoneToZones()[shadPatchI].slaveToMaster
                        (
                            slaveZoneTraction
                        )()
                    )
                );

                // Accumulate the traction on the master patch
                traction() += tractionForThisShadow;

                // Update contactPerShadow field
                // Note: this is used by thermalContact to know which faces
                // are in contact
                const scalarField magTraction(mag(tractionForThisShadow));
                const scalar tol = 1e-6*gMax(magTraction);
                scalarField& contactForThisShadow =
                    contactPerShadow()[shadPatchI];
                forAll(contactForThisShadow, faceI)
                {
                    if (magTraction[faceI] > tol)
                    {
                        contactForThisShadow[faceI] = 1.0;
                    }
                    else
                    {
                        contactForThisShadow[faceI] = 0.0;
                    }
                }
            }
        }
    }
    else
    {
        // Set the traction on the slave patch
        // The master stores the friction and normal models, so we need to find
        // which models correspond to the current shadow
        traction() =
            frictionModelForThisSlave().slaveTraction()
          + normalModelForThisSlave().slavePressure();

        // TESTING - START
        // Scale traction vectors on faces, which share an edge with the
        // downstream patch
        // This is an attempt to fix an issue where the first row of faces
        // deform unphysically when being drawn into the die
        if (scaleFaceTractionsNearDownstreamPatch_)
        {
            traction() *= scaleTractionField();
        }
        // TESTING - END

        // Update contactPerShadow field
        // Note: this is used by thermalContact to know which faces
        // are in contact
        const scalarField magTraction(mag(traction()));
        const scalar tol = 1e-6*gMax(magTraction);
        scalarField& contactForThisShadow = contactPerShadow()[0];
        forAll(contactForThisShadow, faceI)
        {
            if (magTraction[faceI] > tol)
            {
                contactForThisShadow[faceI] = 1.0;
            }
            else
            {
                contactForThisShadow[faceI] = 0.0;
            }
        }
    }

    // Accumulate the contact indicator field
    contact_ = 0.0;
    PtrList<scalarField>& contactPerShadow = this->contactPerShadow();
    forAll(contactPerShadow, shadI)
    {
        contact_ += contactPerShadow[shadI];
    }

    // Scale any face in contact with more than one shadow
    if (gMax(contact_) > (1.0 + SMALL))
    {
        forAll(contact_, faceI)
        {
            if (contact_[faceI] > (1.0 + SMALL))
            {
                // Update the contact weights corresponding to each shadow
                scalar sumContact = 0.0;
                forAll(contactPerShadow, shadI)
                {
                    contactPerShadow[shadI][faceI] /= contact_[faceI];
                    sumContact += contactPerShadow[shadI][faceI];
                }

                if (sumContact > (1.0 + SMALL))
                {
                    FatalErrorInFunction
                        << "There is a problem normalising the contact field"
                        << ", sumContact is: " << sumContact
                        << abort(FatalError);
                }

                // Reset accumulated contact face value to 1.0
                contact_[faceI] = 1.0;
            }
        }
    }

    solidTractionFvPatchVectorField::updateCoeffs();
}


Foam::tmp<Foam::scalarField>
Foam::solidContactFvPatchVectorField::frictionHeatRate() const
{
    // Consider storing frictionHeatRate instead of recalculating multiple times

    if (!master_)
    {
        FatalErrorInFunction
            << "Only master can call frictionHeatRate function!"
            << abort(FatalError);
    }

    // For now, we assume traction is constant over time-step
    // Todo: we should use trapezoidal rule
    vectorField curTraction(patch().size(), vector::zero);

    tmp<scalarField> tfrictionHeatRate
    (
        new scalarField(curTraction.size(), 0.0)
    );
    scalarField& frictionHeatRate = tfrictionHeatRate.ref();

    forAll(shadowPatchNames(), shadPatchI)
    {
        // Calculate slip

        const vectorField slavePatchSlip(frictionModels()[shadPatchI].slip());

        const vectorField slaveZoneSlip
        (
            shadowZones()[shadPatchI].patchFaceToGlobal
            (
                slavePatchSlip
            )
        );

        // Interpolate from slave to master

        const vectorField masterZoneSlip
        (
            zoneToZones()[shadPatchI].slaveToMaster(slaveZoneSlip)
        );

        const vectorField masterPatchSlip
        (
            zone().globalFaceToPatch(masterZoneSlip)
        );

        const scalar deltaT =
            patch().boundaryMesh().mesh().time().deltaTValue();

        // Accumulate frictionHeatRate for each shadow patch

        // Rate of dissipated frictional energy for this timestep
        // The dot product of the traction vectors and the slip vectors gives
        // the dissipated frictional energy per unit area; which is always
        // positive
        frictionHeatRate += mag(traction() & (masterPatchSlip/deltaT));
    }

    return tfrictionHeatRate;
}


Foam::PtrList<Foam::scalarField>&
Foam::solidContactFvPatchVectorField::contactPerShadow()
{
    if (contactPerShadow_.size() == 0)
    {
        calcContactPerShadow();
    }

    return contactPerShadow_;
}


const Foam::PtrList<Foam::scalarField>&
Foam::solidContactFvPatchVectorField::contactPerShadow() const
{
    if (contactPerShadow_.size() == 0)
    {
        calcContactPerShadow();
    }

    return contactPerShadow_;
}


void Foam::solidContactFvPatchVectorField::write(Ostream& os) const
{
    // If the shadowPatchIndices pointer is not set then we will assume that the
    // contact models were not created and nothing has changed; so we will just
    // output the input dict unchanged
    if (shadowPatchNames_.size() == 0)
    {
        // Overwrite fields in the dict
        dictionary& dict = const_cast<dictionary&>(dict_);

        dict.remove("gradient");
        dict.remove("value");
        dict.remove("traction");
        dict.remove("pressure");

        //dict.add("gradient", gradient());
        const vectorField& patchValue = *this;

        // Write the dictionary
        dict_.write(os, false);

        writeEntry(os, "gradient", gradient());
        writeEntry(os, "value", patchValue);
        writeEntry(os, "traction", traction());
        writeEntry(os, "pressure", pressure());

        return;
    }

    solidTractionFvPatchVectorField::write(os);

    writeEntry(os, "master", master_);
    const wordList& shadPatchNames = shadowPatchNames();
    if (shadPatchNames.size() == 1)
    {
        writeEntry(os, "shadowPatch", shadPatchNames[0]);
    }
    else
    {
        writeEntry(os, "shadowPatches", shadPatchNames);
    }

    os.writeKeyword("regionOfInterest")
        << regionOfInterest_ << token::END_STATEMENT << nl;
//     writeEntry(os, "regionOfInterest", regionOfInterest_);
    writeEntry(os, "regionOfInterestTopCorner", regionOfInterestTopCorner_);
    writeEntry(os, "regionOfInterestBottomCorner", regionOfInterestBottomCorner_);
    writeEntry(os, "writeZoneVTK", writeZoneVTK_);
    writeEntry(os, "writePointDistanceFields", writePointDistanceFields_);
    writeEntry
    (
        os,
        "scaleFaceTractionsNearDownstreamPatch",
        scaleFaceTractionsNearDownstreamPatch_
    );

    if (scaleFaceTractionsNearDownstreamPatch_)
    {
        writeEntry
        (
            os,
            "downstreamScaleFactor",
            dict_.lookup<scalar>("downstreamScaleFactor")
        );
        writeEntry
        (
            os,
            "downstreamPatchName",
            dict_.lookup<word>("downstreamPatchName")
        );
    }

    if (master_)
    {
        writeEntry(os, "rigidMaster", rigidMaster_);

        if (shadowPatchNames_.size() == 1)
        {
            writeEntry(os, "normalContactModel", normalModels()[0].type());
            normalModels()[0].writeDict(os);

            writeEntry
            (
                os,
                "frictionContactModel",
                frictionModels()[0].type()
            );
            frictionModels()[0].writeDict(os);

            writeEntry
            (
                os,
                "useNewPointDistanceMethod",
                dict_.lookupOrDefault<Switch>
                (
                    "useNewPointDistanceMethod", false
                )
            );

            writeEntry
            (
                os,
                "projectPointsToPatchBoundary",
                dict_.lookupOrDefault<Switch>
                (
                    "projectPointsToPatchBoundary", false
                )
            );

            writeEntry
            (
                os,
                "checkPointDistanceOrientations",
                dict_.lookupOrDefault<Switch>
                (
                    "checkPointDistanceOrientations", false
                )
            );

            writeEntry
            (
                os,
                "usePrevCandidateMasterNeighbors",
                dict_.lookupOrDefault<Switch>
                (
                    "usePrevCandidateMasterNeighbors", false
                )
            );
        }
        else
        {
            forAll(shadowPatchNames_, shadPatchI)
            {
                os  << patch().name() << "_to_"
                    << shadowPatchNames_[shadPatchI] << "_dict" << nl
                    << '{' << endl;

                writeEntry
                (
                    os,
                    "normalContactModel",
                    normalModels()[shadPatchI].type()
                );
                normalModels()[shadPatchI].writeDict(os);

                writeEntry
                (
                    os,
                    "frictionContactModel",
                    frictionModels()[shadPatchI].type()
                );
                frictionModels()[shadPatchI].writeDict(os);

                os  << '}' << endl;
            }
        }
    }

    if (writeZoneVTK_)
    {
        if
        (
            internalField().name() == "D"
         || internalField().name() == "DD"
        )
        {
            Info<< "Writing deformed zones to VTK" << endl;
            const word timeName =
                patch().boundaryMesh().mesh().time().timeName();

            zone().globalPatch().writeVTK("zone_" + timeName);

            forAll(shadowZones(), shadI)
            {
                shadowZones()[shadI].globalPatch().writeVTK
                (
                    "shadowZone_" + timeName
                );
            }
        }
    }


    // Write out point distance fields for master and slave
    if (writePointDistanceFields_ && master())
    {
        if (normalModels().size() != 1)
        {
            FatalErrorInFunction
                << "The 'writePointDistanceFields' is currently only "
                << "implemented for one-to-one contact"
                << abort(FatalError);
        }

        // Take a reference to the mesh for convenience
        const polyMesh& mesh = patch().patch().boundaryMesh().mesh();

        // Create the point mesh, which is needed for the point field
        pointMesh pMesh(mesh);

        // Create the point distance fields

        pointScalarField dist
        (
            IOobject
            (
                "pointDistance",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            dimensionedScalar("zero", dimless, 0.0)
        );

        pointVectorField distVecs
        (
            IOobject
            (
                "pointDistanceVectors",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            dimensionedVector("zero", dimless, vector::zero)
        );

        // Transfer the patch point distances into the dist point field
        {
            // Lookup the master point distance to intersection
            const scalarField masterpd
            (
                zone().globalPointToPatch
                (
                    zoneToZones()[0].masterPointDistanceToIntersection()
                )
            );
            const vectorField masterpdVecs
            (
                zone().globalPointToPatch
                (
                    zoneToZones()[0].masterPointDistanceVectorsToIntersection()
                )
            );

            const labelList& masterMeshPoints = patch().patch().meshPoints();

            forAll(masterpd, pI)
            {
                const label pointID = masterMeshPoints[pI];
                dist[pointID] = masterpd[pI];
                distVecs[pointID] = masterpdVecs[pI];
            }
        }

        {
            const scalarField slavepd
            (
                shadowZones()[0].globalPointToPatch
                (
                    zoneToZones()[0].slavePointDistanceToIntersection()
                )
            );
            const vectorField slavepdVecs
            (
                shadowZones()[0].globalPointToPatch
                (
                    zoneToZones()[0].slavePointDistanceVectorsToIntersection()
                )
            );

            const labelList& slaveMeshPoints =
                patch().patch().boundaryMesh()
                [
                    shadowPatchIndices()[0]
                ].meshPoints();

            forAll(slavepd, pI)
            {
                const label pointID = slaveMeshPoints[pI];
                dist[pointID] = slavepd[pI];
                distVecs[pointID] = slavepdVecs[pI];
            }
        }

        // Write the field
        InfoInFunction
            << "Writing point distance fields: " << dist.name()
            << " and " << distVecs.name() << endl;
        dist.write();
        distVecs.write();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        solidContactFvPatchVectorField
    );
}


// ************************************************************************* //
