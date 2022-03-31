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
    solidContactFvPatchVectorField

\*---------------------------------------------------------------------------*/

#include "solidContactFvPatchVectorField.H"
#include "pointFields.H"
#include "globalPolyBoundaryMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void
Foam::solidContactFvPatchVectorField::moveZonesToDeformedConfiguration()
{
    // Only the master moves the zones
    if (!master_)
    {
        return;
    }

    // Method
    // We will interpolate the patch face displacements to the patch vertices
    // and then add these vertex/point displacements to the initial patch
    // points
    // We need to take care in parallel, and also realize that the solidModel
    // might have a moving or stationary mesh

    // Assemble the zone face displacement field to move the zones

    // First we will move the master zone

    vectorField zonePointD(zone().globalPatch().points().size(), vector::zero);

    // For a non-moving mesh, we will move the zones by the total
    // displacement, whereas for a moving mesh (updated Lagrangian), we will
    // move the zones by the displacement increment

    if (usePointD_)
    {
        if (movingMesh())
        {
            const pointVectorField& pointDD =
                db().lookupObject<pointVectorField>("pointDD");

            zonePointD = zone().patchPointToGlobal
            (
                pointDD.boundaryField()[patch().index()].patchInternalField()
            );
        }
        else
        {
            const pointVectorField& pointD =
                db().lookupObject<pointVectorField>("pointD");

            zonePointD = zone().patchPointToGlobal
            (
                pointD.boundaryField()[patch().index()].patchInternalField()
            );
        }
    }
    else
    {
        tmp<vectorField> zoneD;
        if (movingMesh())
        {
            const volVectorField& DD =
                db().lookupObject<volVectorField>("DD");

            zoneD = zone().patchFaceToGlobal
            (
                DD.boundaryField()[patch().index()]
            );
        }
        else
        {
            const volVectorField& D = db().lookupObject<volVectorField>("D");
            zoneD = zone().patchFaceToGlobal
            (
                D.boundaryField()[patch().index()]
            );
        }

        //- Interpolated the face values to point values
        zonePointD = zone().interpolator().faceToPointInterpolate(zoneD);
    }

    // The zone deformed points are the initial position plus the
    // displacement
    const pointField zoneNewPoints
    (
        zone().patchPointToGlobal(patch().patch().localPoints())
      + zonePointD
    );

    // Move to new points
    zone().movePoints(zoneNewPoints);


    // Secondly we will move the shadow zones
    forAll(shadowPatchNames(), shadPatchI)
    {
        vectorField shadowZonePointD
        (
            shadowZones()[shadPatchI].globalPatch().points().size(),
            vector::zero
        );

        if (usePointD_)
        {
            if (movingMesh())
            {
                const pointVectorField& pointDD =
                    db().lookupObject<pointVectorField>("pointDD");

                shadowZonePointD =
                    shadowZones()[shadPatchI].patchPointToGlobal
                    (
                        pointDD.boundaryField()
                        [
                            shadowPatchIndices()[shadPatchI]
                        ].patchInternalField()
                    );
            }
            else
            {
                const pointVectorField& pointD =
                    db().lookupObject<pointVectorField>("pointD");

                shadowZonePointD =
                    shadowZones()[shadPatchI].patchPointToGlobal
                    (
                        pointD.boundaryField()
                        [
                            shadowPatchIndices()[shadPatchI]
                        ].patchInternalField()
                    );
            }
        }
        else
        {
            tmp<vectorField> shadowZoneD;
            if (movingMesh())
            {
                const volVectorField& DD =
                    db().lookupObject<volVectorField>("DD");

                shadowZoneD =
                    shadowZones()[shadPatchI].patchFaceToGlobal
                    (
                        DD.boundaryField()[shadowPatchIndices()[shadPatchI]]
                    );
            }
            else
            {
                const volVectorField& D =
                    db().lookupObject<volVectorField>("D");

                shadowZoneD =
                    shadowZones()[shadPatchI].patchFaceToGlobal
                    (
                        D.boundaryField()[shadowPatchIndices()[shadPatchI]]
                    );
            }

            // Interpolate face values to points
            shadowZonePointD =
                shadowZones()[shadPatchI].interpolator().faceToPointInterpolate
                (
                    shadowZoneD
                );

        }

        // The zone deformed points are the initial position plus the
        // displacement
        const pointField shadowZoneNewPoints
        (
            shadowZones()[shadPatchI].patchPointToGlobal
            (
                shadowZones()[shadPatchI].patch().localPoints()
            )
          + shadowZonePointD
        );

        // Move shadow patch to the new points
        shadowZones()[shadPatchI].movePoints(shadowZoneNewPoints);
    }
}


void Foam::solidContactFvPatchVectorField::calcZone() const
{
    if (debug)
    {
        InfoInFunction
            << patch().name() << " : making the zone" << endl;
    }

    if (!master_)
    {
        FatalErrorInFunction
            << "Trying to create zone on a slave" << abort(FatalError);
    }

    if (zonePtr_)
    {
        FatalErrorInFunction
            << "pointer already set" << abort(FatalError);
    }

    // Note: the main mesh will either be in the initial configuration or the
    // updated configuration
    zonePtr_ = new globalPolyPatch(patch().patch());
}


void Foam::solidContactFvPatchVectorField::calcShadowZones() const
{
    if (debug)
    {
        InfoInFunction
            << patch().name() << " : making the shadow zones" << endl;
    }

    if (!master_)
    {
        FatalErrorInFunction
            << "Trying to create shadow zones on a slave" << abort(FatalError);
    }

    if (!shadowZones_.empty())
    {
        FatalErrorInFunction
            << "pointer already set" << abort(FatalError);
    }

    const wordList& shadPatchNames = shadowPatchNames();

    shadowZones_.setSize(shadPatchNames.size());

    const polyBoundaryMesh& bMesh = patch().patch().boundaryMesh();

    forAll(shadowZones_, shadPatchI)
    {
        // Note: the main mesh will either be in the initial configuration or
        // the updated configuration
        shadowZones_.set
        (
            shadPatchI,
            new globalPolyPatch(bMesh[shadPatchNames[shadPatchI]])
        );
    }
}


void Foam::solidContactFvPatchVectorField::calcZoneToZones() const
{
    if (debug)
    {
        InfoInFunction
            << patch().name() << " : making the zoneToZone" << endl;
    }

    // Create zone-to-zone interpolation
    if (!zoneToZones_.empty())
    {
        FatalErrorInFunction
            << "Zone to zone interpolation already calculated"
            << abort(FatalError);
    }

    // Check master and slave patch
    const volVectorField& field =
        db().lookupObject<volVectorField>
        (
            this->internalField().name()
        );

    zoneToZones_.setSize(shadowPatchNames().size());

    forAll (zoneToZones_, shadPatchI)
    {
        const solidContactFvPatchVectorField& shadowPatchField =
            refCast<const solidContactFvPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndices()[shadPatchI]]
            );

        if (master_)
        {
            if (shadowPatchField.master() == true)
            {
                FatalErrorInFunction
                    << "There are two master patches!" << abort(FatalError);
            }
        }
        else
        {
            if (shadowPatchField.master() == false)
            {
                FatalErrorInFunction
                    << "There is no master patch!" << abort(FatalError);
            }
        }

        if (master_)
        {
            // Create interpolation for patches
            zoneToZones_.set
            (
                shadPatchI,
                new ggiStandAlonePatchInterpolation
                (
                    zone().globalPatch(),
                    shadowZones()[shadPatchI].globalPatch(),
                    tensorField(0),
                    tensorField(0),
                    vectorField(0), // Slave-to-master separation. Bug fix
                    true,           // global data
                    0,              // Master non-overlapping face tolerances
                    0,              // Slave non-overlapping face tolerances
                    // Do not rescale weighting factors, as it is wrong on
                    // partially covered faces
                    false,
                    quickReject_,
                    regionOfInterest_
                )
            );

            // Check which point distance calculation method to use
            const Switch useNewPointDistanceMethod =
                dict_.lookupOrDefault<Switch>
                (
                    "useNewPointDistanceMethod", false
                );

            Info<< "    " << type() << ": " << patch().name() << nl
                << "        useNewPointDistanceMethod: "
                << useNewPointDistanceMethod
                << endl;

            zoneToZones_[shadPatchI].useNewPointDistanceMethod() =
                useNewPointDistanceMethod;

            // Check if the projectPointsToPatchBoundary switch is set
            const Switch projectPointsToPatchBoundary =
                dict_.lookupOrDefault<Switch>
                (
                    "projectPointsToPatchBoundary",
                    false
                );

            Info<< "        projectPointsToPatchBoundary: "
                << projectPointsToPatchBoundary
                << endl;

            zoneToZones_[shadPatchI].projectPointsToPatchBoundary() =
                projectPointsToPatchBoundary;

            if (dict_.found("checkPointDistanceOrientations"))
            {
                const Switch checkPointDistanceOrientations =
                    Switch(dict_.lookup("checkPointDistanceOrientations"));

                Info<< "        checkPointDistanceOrientations: "
                    << checkPointDistanceOrientations
                    << endl;

                zoneToZones_[shadPatchI].checkPointDistanceOrientations() =
                    checkPointDistanceOrientations;
            }

            // Check if the usePrevCandidateMasterNeighbors switch is set
            const Switch usePrevCandidateMasterNeighbors =
                dict_.lookupOrDefault<Switch>
                (
                    "usePrevCandidateMasterNeighbors",
                    false
                );

            Info<< "        usePrevCandidateMasterNeighbors: "
                << usePrevCandidateMasterNeighbors
                << endl;

            zoneToZones_[shadPatchI].usePrevCandidateMasterNeighbors() =
                usePrevCandidateMasterNeighbors;
        }
        else
        {
            FatalErrorInFunction
                << "Attempting to create GGIInterpolation on a shadow patch"
                << abort(FatalError);
        }
    }
}


void Foam::solidContactFvPatchVectorField::calcContactPerShadow() const
{
    if (contactPerShadow_.size() > 0)
    {
        FatalErrorInFunction
            << "already calculated"
            << abort(FatalError);
    }

    contactPerShadow_.setSize(shadowPatchNames().size());

    forAll(contactPerShadow_, i)
    {
        contactPerShadow_.set
        (
            i,
            new scalarField(patch().size(), 0.0)
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::globalPolyPatch&
Foam::solidContactFvPatchVectorField::zone() const
{
    const globalPolyBoundaryMesh& gpbMesh
    (
        globalPolyBoundaryMesh::New(patch().boundaryMesh().mesh())
    );
    if (master_)
    {
        return gpbMesh[patch().patch()];
    }
    else
    {
        return gpbMesh[shadowPatchNames()[0]];
    }

//     if (master_)
//     {
//         if (!zonePtr_)
//         {
//             calcZone();
//         }
//
//         return *zonePtr_;
//     }
//     else
//     {
//         const volVectorField& field =
//             db().lookupObject<volVectorField>
//             (
//                 this->internalField().name()
//             );
//
//         const solidContactFvPatchVectorField& shadowPatchField =
//             refCast<const solidContactFvPatchVectorField>
//             (
//                 field.boundaryField()[shadowPatchIndices()[0]]
//             );
//
//         return shadowPatchField.zone();
//     }
}


Foam::globalPolyPatch& Foam::solidContactFvPatchVectorField::zone()
{
    const globalPolyBoundaryMesh& gpbMesh
    (
        globalPolyBoundaryMesh::New(patch().boundaryMesh().mesh())
    );
    if (master_)
    {
        return const_cast<globalPolyPatch&>(gpbMesh[patch().patch()]);
    }
    else
    {
        return const_cast<globalPolyPatch&>(gpbMesh[shadowPatchNames()[0]]);
    }
//     if (master_)
//     {
//         if (!zonePtr_)
//         {
//             calcZone();
//         }
//
//         return *zonePtr_;
//     }
//     else
//     {
//         const volVectorField& field =
//             db().lookupObject<volVectorField>
//             (
//                 this->internalField().name()
//             );
//
//         solidContactFvPatchVectorField& shadowPatchField =
//             const_cast<solidContactFvPatchVectorField&>
//             (
//                 refCast<const solidContactFvPatchVectorField>
//                 (
//                     field.boundaryField()[shadowPatchIndices()[0]]
//                 )
//             );
//
//         return shadowPatchField.zone();
//     }
}


const Foam::PtrList<Foam::globalPolyPatch>&
Foam::solidContactFvPatchVectorField::shadowZones() const
{
    if (master_)
    {
        if (shadowZones_.empty())
        {
            calcShadowZones();
        }

        return shadowZones_;
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

        return shadowPatchField.shadowZones();
    }
}


Foam::PtrList<Foam::globalPolyPatch>&
Foam::solidContactFvPatchVectorField::shadowZones()
{
    if (master_)
    {
        if (shadowZones_.empty())
        {
            calcShadowZones();
        }

        return shadowZones_;
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

        return shadowPatchField.shadowZones();
    }
}


const Foam::PtrList<Foam::ggiStandAlonePatchInterpolation>&
Foam::solidContactFvPatchVectorField::zoneToZones() const
{
    if (master_)
    {
        if (zoneToZones_.empty())
        {
            calcZoneToZones();
        }

        return zoneToZones_;
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

        return shadowPatchField.zoneToZones();
    }
}


Foam::PtrList<Foam::ggiStandAlonePatchInterpolation>&
Foam::solidContactFvPatchVectorField::zoneToZones()
{
    if (master_)
    {
        if (zoneToZones_.empty())
        {
            calcZoneToZones();
        }

        return zoneToZones_;
    }
    else
    {
        // We will const_cast the shadow patch so we can delete the weights when
        // the zones move
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

        return shadowPatchField.zoneToZones();
    }

}


const Foam::ggiStandAlonePatchInterpolation&
Foam::solidContactFvPatchVectorField::zoneToZoneForThisSlave() const
{
    if (master_)
    {
        FatalErrorInFunction
            << "The master patch is not allowed to call this function"
            << abort(FatalError);
    }

    // The master may have multiple slaves so we need to find which zoneToZone
    // corresponds to the current slave patch
    const wordList& shadPatchNames = shadowPatchField().shadowPatchNames();
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

    // Return the zoneToZone between the master and the current patch
    return zoneToZones()[masterShadowID];
}


const Foam::globalPolyPatch&
Foam::solidContactFvPatchVectorField::zoneForThisSlave() const
{
    if (master_)
    {
        FatalErrorInFunction
            << "The master patch is not allowed to call this function"
            << abort(FatalError);
    }

    // The master may have multiple slaves so we need to find which shadowZone
    // corresponds to the current slave patch
    const wordList& shadPatchNames = shadowPatchField().shadowPatchNames();
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

    // Return the zoneToZone between the master and the current patch
    return shadowZones()[masterShadowID];
}


// ************************************************************************* //
