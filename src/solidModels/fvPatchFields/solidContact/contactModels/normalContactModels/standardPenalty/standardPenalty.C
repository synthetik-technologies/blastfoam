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
    standardPenalty

\*---------------------------------------------------------------------------*/

#include "standardPenalty.H"
#include "volFields.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

  defineTypeNameAndDebug(standardPenalty, 0);
  addToRunTimeSelectionTable(normalContactModel, standardPenalty, dictionary);


// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //


void standardPenalty::calcPenaltyFactor() const
{
    if (penaltyFactor_ > -SMALL)
    {
        FatalErrorInFunction
            << "value already set!" << abort(FatalError);
    }

    // Set penalty factor

    // Approximate penaltyFactor from the mechanical properties
    // This can then be scaled using the penaltyScale

    const label masterPatchIndex =  masterPatchID();
    const label slavePatchIndex =  slavePatchID();

    // Lookup implicit stiffness = 2*mu + lambda, approximately equal to the
    // bulk modulus
    const volScalarField& impK = mesh_.lookupObject<volScalarField>("impK");

    // Note: for solidRigidContact, only the master index is set
    if (masterPatchIndex > -1 && slavePatchIndex > -1)
    {
        // Avarage contact patch bulk modulus
        const scalar masterK = gAverage(impK.boundaryField()[masterPatchIndex]);
        const scalar slaveK = gAverage(impK.boundaryField()[slavePatchIndex]);

        // Simple average
        //const scalar bulkModulus = 0.5*(masterK + slaveK);
        // Harmonic average
        const scalar bulkModulus = 2.0*(masterK*slaveK)/(masterK + slaveK);

        // Average contact patch face area
        const scalar masterMagSf =
            gAverage(mesh_.magSf().boundaryField()[masterPatchIndex]);
        const scalar slaveMagSf =
            gAverage(mesh_.magSf().boundaryField()[slavePatchIndex]);
        const scalar faceArea = 0.5*(masterMagSf + slaveMagSf);

        // Average contact patch cell volume
        scalarField masterV(mesh_.boundary()[masterPatchIndex].size(), 0.0);
        scalarField slaveV(mesh_.boundary()[slavePatchIndex].size(), 0.0);
        const volScalarField::Internal& V = mesh_.V();
        {
            const unallocLabelList& faceCells =
                mesh_.boundary()[masterPatchIndex].faceCells();
            forAll(mesh_.boundary()[masterPatchIndex], facei)
            {
                masterV[facei] = V[faceCells[facei]];
            }
        }
        {
            const unallocLabelList& faceCells =
                mesh_.boundary()[slavePatchIndex].faceCells();
            forAll(mesh_.boundary()[slavePatchIndex], facei)
            {
                slaveV[facei] = V[faceCells[facei]];
            }
        }
        const scalar cellVolume = 0.5*(gAverage(masterV) + gAverage(slaveV));

        // Approximate penalty factor based on Hallquist et al.
        // we approximate penalty factor for traction instead of force
        penaltyFactor_ = penaltyScale_*bulkModulus*faceArea/cellVolume;
    }
    else if (slavePatchIndex > -1)
    {
        // Avarage contact patch bulk modulus
        const scalar bulkModulus =
            gAverage(impK.boundaryField()[slavePatchIndex]);

        // Average contact patch face area
        const scalar faceArea =
            gAverage(mesh_.magSf().boundaryField()[slavePatchIndex]);

        // Average contact patch cell volume
        scalarField slaveV(mesh_.boundary()[slavePatchIndex].size(), 0.0);
        const volScalarField::Internal& V = mesh_.V();
        {
            const unallocLabelList& faceCells =
                mesh_.boundary()[slavePatchIndex].faceCells();
            forAll(mesh_.boundary()[slavePatchIndex], facei)
            {
                slaveV[facei] = V[faceCells[facei]];
            }
        }
        const scalar cellVolume = gAverage(slaveV);

        // Approximate penalty factor based on Hallquist et al.
        // we approximate penalty factor for traction instead of force
        penaltyFactor_ = penaltyScale_*bulkModulus*faceArea/cellVolume;
    }
    else
    {
        FatalErrorIn("void standardPenalty::calcPenaltyFactor() const")
            << "This is unexpected! Neither the master nor slave index are set!"
            << abort(FatalError);
    }

    Info<< "    normal penalty factor: " << penaltyFactor_ << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

standardPenalty::standardPenalty
(
    const word& name,
    const fvPatch& patch,
    const dictionary& dict,
    const label masterPatchID,
    const label slavePatchID,
    const standAlonePatch& masterFaceZonePatch,
    const standAlonePatch& slaveFaceZonePatch
)
:
    normalContactModel
    (
        name,
        patch,
        dict,
        masterPatchID,
        slavePatchID,
        masterFaceZonePatch,
        slaveFaceZonePatch
    ),
    normalContactModelDict_(dict.subDict(name + "NormalModelDict")),
    mesh_(patch.boundaryMesh().mesh()),
    slavePressureVolField_
    (
        IOobject
        (
            "slavePressure_" + mesh_.boundaryMesh()[slavePatchID].name(),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", dimPressure, vector::zero)
    ),
    areaInContactVolField_
    (
        IOobject
        (
            "areaInContact_" + mesh_.boundaryMesh()[slavePatchID].name(),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimArea, 0.0)
    ),
    penaltyFactor_(-1),
    penaltyScale_
    (
        normalContactModelDict_.lookupOrDefault<scalar>("penaltyScale", 1.0)
    ),
    relaxFac_
    (
        normalContactModelDict_.lookupOrDefault<scalar>
        (
            "relaxationFactor", 0.02
        )
    ),
    averagePenetration_(0),
    minPenetration_(0),
    epsilon0_
    (
        normalContactModelDict_.lookupOrDefault<scalar>("epsilon0", 0.0)
    ),
    contactIterNum_(0)
{}


standardPenalty::standardPenalty(const standardPenalty& nm)
:
    normalContactModel(nm),
    normalContactModelDict_(nm.normalContactModelDict_),
    mesh_(nm.mesh_),
    slavePressureVolField_(nm.slavePressureVolField_),
    areaInContactVolField_(nm.areaInContactVolField_),
    penaltyFactor_(nm.penaltyFactor_),
    penaltyScale_(nm.penaltyScale_),
    relaxFac_(nm.relaxFac_),
    averagePenetration_(nm.averagePenetration_),
    minPenetration_(nm.minPenetration_),
    epsilon0_(nm.epsilon0_),
    contactIterNum_(nm.contactIterNum_)
{}


// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //


void standardPenalty::correct
(
    const vectorField& slavePatchFaceNormals,
    const scalarField& slavePointPenetration,
    const vectorField& slaveDU,
    const vectorField& masterDUInterpToSlave
)
{
    // Preliminaries
    const fvMesh& mesh = mesh_;
    const label slavePatchIndex = slavePatchID();

    // Calculate area in contact for slave patch

    const faceList& slavePatchLocalFaces =
        mesh.boundaryMesh()[slavePatchIndex].localFaces();

    const pointField& slavePatchLocalPoints =
        mesh.boundaryMesh()[slavePatchIndex].localPoints();

    scalarField slavePatchLocalFaceAreas(slavePatchLocalFaces.size(), 0.0);

    scalarField& areaInContact = this->areaInContact();
    forAll(slavePatchLocalFaces, faceI)
    {
        areaInContact[faceI] =
            slavePatchLocalFaces[faceI].areaInContact
            (
                slavePatchLocalPoints,
                slavePointPenetration
            );

        if (areaInContact[faceI] < -SMALL)
        {
            const labelList& labels = slavePatchLocalFaces[faceI];
            scalarField vertexValue(labels.size());
            forAll(labels, i)
            {
                vertexValue[i] = slavePointPenetration[labels[i]];
            }

            FatalErrorIn(type())
                << "areaInContact is less than zero!" << nl
                << "areaInContact[" << faceI << "] = " << areaInContact[faceI]
                << nl
                << "vertexValue = " << vertexValue << nl
                << endl;
        }

        slavePatchLocalFaceAreas[faceI] =
            mag(slavePatchLocalFaces[faceI].normal(slavePatchLocalPoints));
    }

    // Calculate the point pressures
    // We will also record the average and minium penetrations

    const scalar penaltyFac = penaltyFactor();
    scalarField totalSlavePointPressure(slavePointPenetration.size(), 0.0);
    averagePenetration_ = 0.0;
    minPenetration_ = 0.0;
    int nPointsInContact = 0;

    forAll(totalSlavePointPressure, pointI)
    {
        // Take copy of penetration
        const scalar d = slavePointPenetration[pointI];

        // Note: penetration is negative for points in contact
        {
            // The force is linearly proportional the penetration, like a spring
            if (d < epsilon0_)
            {
                totalSlavePointPressure[pointI] =
                    max(penaltyFac*(epsilon0_ - d), 0.0);

                averagePenetration_ += d;
                minPenetration_ = min(minPenetration_, d);
                nPointsInContact++;
            }
            else
            {
                totalSlavePointPressure[pointI] = 0.0;
            }
        }
    }

    // Find global minimum penetration
    // IB 11/2018
    reduce(minPenetration_, minOp<scalar>());

    // Update the average penetration
    reduce(averagePenetration_, sumOp<scalar>());
    reduce(nPointsInContact, sumOp<label>());
    if (nPointsInContact > 0)
    {
        averagePenetration_ /= nPointsInContact;
    }
    else
    {
        averagePenetration_ = 0.0;
    }


    // Interpolate point pressures to the faces

    // Create local patch interpolation: No need to interpolate using the entire
    // face zone patch
    primitivePatchInterpolation localSlaveInterpolator
    (
        mesh.boundaryMesh()[slavePatchIndex]
    );

    // Interpolate point pressures to the face centres and apply in the negative
    // normal direction
    vectorField newSlaveTraction
    (
        localSlaveInterpolator.pointToFaceInterpolate<scalar>
        (
            totalSlavePointPressure
        )*(-slavePatchFaceNormals)
    );

    // Under-relax pressure/traction
    // Note: slavePressure_ is really a traction vector
    slavePressure() =
        relaxFac_*newSlaveTraction + (1.0 - relaxFac_)*slavePressure();
}


scalar standardPenalty::penaltyFactor() const
{
    if (penaltyFactor_ < -SMALL)
    {
        calcPenaltyFactor();
    }

    return penaltyFactor_;
}


scalar standardPenalty::updatePenaltyScale(const scalar previousPenaltyScale)
{
    if (previousPenaltyScale > 0.0)
    {
        // Lookup initial value for penaltyScale
        const scalar initialPenaltyScale =
            normalContactModelDict_.lookupOrDefault<scalar>
            (
                "penaltyScale", 1.0
            );

        if (mag(initialPenaltyScale - penaltyScale_) < SMALL)
        {
            // After a topo change, use the previous value of penalty scale
            penaltyScale_ = previousPenaltyScale;
        }
    }

    return penaltyScale_;
}


void standardPenalty::autoMap(const fvPatchFieldMapper& m)
{
    if (debug)
    {
        InfoIn
        (
            "void standardPenalty::autoMap(const fvPatchFieldMapper& m)"
        )   << "autoMap" << endl;
    }

    normalContactModel::autoMap(m);

    // The internal fields for the volFields should always be zero
    // We will reset them as they may not be zero after field advection
    slavePressureVolField_.primitiveFieldRef() = vector::zero;
    areaInContactVolField_.primitiveFieldRef() = 0.0;
}


void standardPenalty::writeDict(Ostream& os) const
{
    // Update the penalty scale in the dictionary
    normalContactModelDict_.set("penaltyScale", penaltyScale_);
    writeEntry(os, name() + "NormalModelDict", normalContactModelDict_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
