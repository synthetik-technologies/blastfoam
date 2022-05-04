/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dynMeshTools.H"
#include "polyMesh.H"
#include "hexMatcher.H"
#include "faceZone.H"
#include "syncTools.H"
#include "polyTopoChange.H"
#include "polyAddFace.H"
#include "polyModifyFace.H"
#include "removeCells.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "fvMeshTools.H"

#include "polyModifyFace.H"
#include "polyAddFace.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meshTools::getFaceInfo
(
    const polyMesh& mesh,
    const label faceI,
    label& patchID,
    label& zoneID,
    label& zoneFlip
)
{
    patchID = -1;

    if (!mesh.isInternalFace(faceI))
    {
        patchID = mesh.boundaryMesh().whichPatch(faceI);
    }

    zoneID = mesh.faceZones().whichZone(faceI);

    zoneFlip = false;

    if (zoneID > -1)
    {
        const faceZone& fZone = mesh.faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
    }
}


Foam::label Foam::meshTools::addFace
(
    polyTopoChange& meshMod,
    const polyMesh& mesh,
    const label faceI,
    const face& newFace,
    const label own,
    const label nei
)
{
    // Set face information
    label patchID, zoneID, zoneFlip;
    meshTools::getFaceInfo(mesh, faceI, patchID, zoneID, zoneFlip);

    // Set new face index to -1
    label newFaceI = -1;

    if ((nei == -1) || (own < nei))
    {
        // Ordering is ok, add the face
        newFaceI = meshMod.setAction
        (
            polyAddFace
            (
                newFace,                    // face
                own,                        // owner
                nei,                        // neighbour
                -1,                         // master point
                -1,                         // master edge
                faceI,                      // master face for addition
                false,                      // flux flip
                patchID,                    // patch for face
                zoneID,                     // zone for face
                zoneFlip                    // face zone flip
            )
        );
    }
    else
    {
        // Ordering is flipped, reverse face and flip owner/neighbour
        newFaceI = meshMod.setAction
        (
            polyAddFace
            (
                newFace.reverseFace(),      // face
                nei,                        // owner
                own,                        // neighbour
                -1,                         // master point
                -1,                         // master edge
                faceI,                      // master face for addition
                false,                      // flux flip
                patchID,                    // patch for face
                zoneID,                     // zone for face
                zoneFlip                    // face zone flip
            )
        );
    }

    return newFaceI;
}


Foam::label Foam::meshTools::addInternalFace
(
    polyTopoChange& meshMod,
    const polyMesh& mesh,
    const label meshFaceI,
    const label meshPointI,
    const face& newFace,
    const label own,
    const label nei
)
{
    // Check whether this is an internal face
    if (mesh.isInternalFace(meshFaceI))
    {
        return meshMod.setAction
        (
            polyAddFace
            (
                newFace,                    // face
                own,                        // owner
                nei,                        // neighbour
                -1,                         // master point
                -1,                         // master edge
                meshFaceI,                  // master face for addition
                false,                      // flux flip
                -1,                         // patch for face
                -1,                         // zone for face
                false                       // face zone flip
            )
        );
    }
    else
    {
        // This is not an internal face. Add face out of nothing
        return meshMod.setAction
        (
            polyAddFace
            (
                newFace,                    // face
                own,                        // owner
                nei,                        // neighbour
                -1,                         // master point
                -1,                         // master edge
                -1,                         // master face for addition
                false,                      // flux flip
                -1,                         // patch for face
                -1,                         // zone for face
                false                       // face zone flip
            )
        );
    }
}


void Foam::meshTools::modifyFace
(
    polyTopoChange& meshMod,
    const polyMesh& mesh,
    const label faceI,
    const face& newFace,
    const label own,
    const label nei
)
{
    // Set face inforomation
    label patchID, zoneID, zoneFlip;
    meshTools::getFaceInfo(mesh, faceI, patchID, zoneID, zoneFlip);

    // Get owner/neighbour addressing and mesh faces
    const labelList& owner = mesh.faceOwner();
    const labelList& neighbour = mesh.faceNeighbour();

    const faceList& meshFaces = mesh.faces();

    if
    (
        (own != owner[faceI])
     || (
            mesh.isInternalFace(faceI)
         && (nei != neighbour[faceI])
        )
     || (newFace != meshFaces[faceI])
    )
    {
        // Either:
        // 1. Owner index does not correspond to mesh owner,
        // 2. Neighbour index does not correspond to mesh neighbour,
        // 3. New face does not correspond to mesh face
        // So we need to modify this face
        if ((nei == -1) || (own < nei))
        {
            // Ordering is ok, add the face
            meshMod.setAction
            (
                polyModifyFace
                (
                    newFace,            // modified face
                    faceI,              // label of face being modified
                    own,                // owner
                    nei,                // neighbour
                    false,              // face flip
                    patchID,            // patch for face
                    false,              // remove from zone
                    zoneID,             // zone for face
                    zoneFlip            // face flip in zone
                )
            );
        }
        else
        {
            // Ordering is flipped, reverse face and flip owner/neighbour
            meshMod.setAction
            (
                polyModifyFace
                (
                    newFace.reverseFace(),  // modified face
                    faceI,                  // label of face being modified
                    nei,                    // owner
                    own,                    // neighbour
                    false,                  // face flip
                    patchID,                // patch for face
                    false,                  // remove from zone
                    zoneID,                 // zone for face
                    zoneFlip                // face flip in zone
                )
            );
        }
    }
}


void Foam::meshTools::changePatchFace
(
    polyTopoChange& meshMod,
    const polyMesh& mesh,
    const label faceI,
    const label newPatchID
)
{
    // Set face inforomation
    label patchID, zoneID, zoneFlip;
    meshTools::getFaceInfo(mesh, faceI, patchID, zoneID, zoneFlip);

    meshMod.setAction
    (
        polyModifyFace
        (
            mesh.faces()[faceI],    // modified face
            faceI,                  // label of face being modified
            mesh.faceOwner()[faceI],// owner
            -1,                     // neighbour
            false,                  // face flip
            newPatchID,             // patch for face
            false,                  // remove from zone
            zoneID,                 // zone for face
            zoneFlip                // face flip in zone
        )
    );
}


void Foam::meshTools::checkInternalOrientation
(
    const polyTopoChange& meshMod,
    const label cellI,
    const label faceI,
    const point& ownPt,
    const point& neiPt,
    const face& newFace
)
{
    const face compactFace(identity(newFace.size()));

    // Get compact points
    const pointField compactPoints(meshMod.points(), newFace);

    const vector n(compactFace.normal(compactPoints));
    const vector dir(neiPt - ownPt);

    // Check orientation error
    if ((dir & n) < 0)
    {
        FatalErrorInFunction
            << "cell:" << cellI << nl
            << "face:" << faceI << nl
            << "newFace:" << newFace << nl
            << "coords:" << compactPoints << nl
            << "ownPt:" << ownPt << nl
            << "neiPt:" << neiPt << nl
            << abort(FatalError);
    }

    // Note: report significant non-orthogonality error
    const scalar severeNonOrthogonalityThreshold =
        cos(70.0/180.0*constant::mathematical::pi);

    const vector fcToOwn(compactFace.centre(compactPoints) - ownPt);

    const scalar s = (fcToOwn & n)/(mag(fcToOwn) + VSMALL);

    if (s > severeNonOrthogonalityThreshold)
    {
        WarningInFunction
            << "Detected severely non-orthogonal face with non-orthogonality: "
            << acos(s)/constant::mathematical::pi*180.0 << nl
            << "cell:" << cellI << " old face:" << faceI << nl
            << "newFace: " << newFace << nl
            << "coords: " << compactPoints << nl
            << "ownPt: " << ownPt << nl
            << "neiPt: " << neiPt << nl
            << "s: " << s << nl
            << endl;
    }
}


void Foam::meshTools::checkBoundaryOrientation
(
    const polyTopoChange& meshMod,
    const label cellI,
    const label faceI,
    const point& ownPt,
    const point& boundaryPt,
    const face& newFace
)
{
    face compactFace(identity(newFace.size()));
    pointField compactPoints(meshMod.points(), newFace);

    vector n(compactFace.normal(compactPoints));

    vector dir(boundaryPt - ownPt);

    if ((dir & n) < 0)
    {
        FatalErrorInFunction
            << "cell:" << cellI << " old face:" << faceI
            << " newFace:" << newFace
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " boundaryPt:" << boundaryPt
            << abort(FatalError);
    }

    // Note: report significant non-orthogonality error
    const scalar severeNonOrthogonalityThreshold =
        cos(70.0/180.0*constant::mathematical::pi);

    const vector fcToOwn(compactFace.centre(compactPoints) - ownPt);

    const scalar s = (fcToOwn & n)/(mag(fcToOwn) + VSMALL);

    if (s > severeNonOrthogonalityThreshold)
    {
        WarningInFunction
            << "Detected severely non-orthogonal face with non-orthogonality: "
            << acos(s)/constant::mathematical::pi*180.0
            << "cell:" << cellI << " old face:" << faceI
            << " newFace:" << newFace
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " boundaryPt:" << boundaryPt
            << " s:" << s
            << endl;
    }
}


void Foam::meshTools::checkFaceOrientation
(
    const polyTopoChange& meshMod,
    const polyMesh& mesh,
    const label& faceI,
    const face& newFace
)
{
    const vectorField& C = mesh.cellCentres();
    if (mesh.isInternalFace(faceI))
    {
        // Get old owner/neighbour indices
        const label own = mesh.faceOwner()[faceI];
        const label nei = mesh.faceNeighbour()[faceI];

        meshTools::checkInternalOrientation
        (
            meshMod,
            own,
            faceI,
            C[own],
            C[nei],
            newFace
        );
    }
    else
    {
        // Get face centres and old owner
        const vectorField& Cf = mesh.faceCentres();
        const label own = mesh.faceOwner()[faceI];

        meshTools::checkBoundaryOrientation
        (
            meshMod,
            own,
            faceI,
            C[own],
            Cf[faceI],
            newFace
        );
    }
}


Foam::label Foam::meshTools::addPatch
(
    polyMesh& mesh,
    const word& patchName,
    const wordList& groupNames,
    const dictionary& patchDict
)
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    if (pbm.findPatchID(patchName) == -1)
    {
        autoPtr<polyPatch> ppPtr
        (
            polyPatch::New
            (
                patchName,
                patchDict,
                0,
                pbm
            )
        );
        polyPatch& pp = ppPtr();

        forAll(groupNames, i)
        {
            if (!pp.inGroup(groupNames[i]))
            {
                pp.inGroups().append(groupNames[i]);
            }
        }


        // Add patch, create calculated everywhere
        mesh.addPatch
        (
            pbm.size(),
            pp,
            dictionary(),   // do not set specialised patchFields
            calculatedFvPatchField<scalar>::typeName,
            true            // parallel sync'ed addition
        );
    }
    else
    {
        Info<< "Patch '" << patchName
            << "' already exists.  Only "
            << "moving patch faces - type will remain the same"
            << endl;
    }

    return pbm.findPatchID(patchName);
}


Foam::label Foam::meshTools::addPatch
(
    polyMesh& mesh,
    const word& patchName,
    const word& groupName,
    const dictionary& patchDict
)
{
    return addPatch(mesh, patchName, wordList(1, groupName), patchDict);
}

void Foam::meshTools::modifyOrAddFace
(
    polyTopoChange& meshMod,
    const face& f,
    const label facei,
    const label own,
    const bool flipFaceFlux,
    const label newPatchi,
    const label zoneID,
    const bool zoneFlip,

    PackedBoolList& modifiedFace
)
{
    if (!modifiedFace.get(facei))
    {
        // First usage of face. Modify.
        meshMod.setAction
        (
            polyModifyFace
            (
                f,                          // modified face
                facei,                      // label of face
                own,                        // owner
                -1,                         // neighbour
                flipFaceFlux,               // face flip
                newPatchi,                  // patch for face
                false,                      // remove from zone
                zoneID,                     // zone for face
                zoneFlip                    // face flip in zone
            )
        );
        modifiedFace.set(facei);
    }
    else
    {
        // Second or more usage of face. Add.
        meshMod.setAction
        (
            polyAddFace
            (
                f,                          // modified face
                own,                        // owner
                -1,                         // neighbour
                -1,                         // master point
                -1,                         // master edge
                facei,                      // master face
                flipFaceFlux,               // face flip
                newPatchi,                  // patch for face
                zoneID,                     // zone for face
                zoneFlip                    // face flip in zone
            )
        );
    }
}


Foam::label Foam::meshTools::createBaffleFaces
(
    const bool internalFacesOnly,
    const polyMesh& mesh,
    const faceZone& fZone,
    const labelList& newMasterPatches,
    const labelList& newSlavePatches,
    polyTopoChange& meshMod,
    PackedBoolList& modifiedFace
)
{
    label nModified = 0;
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    forAll(newMasterPatches, i)
    {
        // Pass 1. Do selected side of zone
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            label zoneFacei = fZone.whichFace(facei);

            if (zoneFacei != -1)
            {
                if (!fZone.flipMap()[zoneFacei])
                {
                    // Use owner side of face
                    modifyOrAddFace
                    (
                        meshMod,
                        mesh.faces()[facei],    // modified face
                        facei,                  // label of face
                        mesh.faceOwner()[facei],// owner
                        false,                  // face flip
                        newMasterPatches[i],    // patch for face
                        fZone.index(),          // zone for face
                        false,                  // face flip in zone
                        modifiedFace            // modify or add status
                    );
                }
                else
                {
                    // Use neighbour side of face.
                    // To keep faceZone pointing out of original neighbour
                    // we don't need to set faceFlip since that cell
                    // now becomes the owner
                    modifyOrAddFace
                    (
                        meshMod,
                        mesh.faces()[facei].reverseFace(),  // modified face
                        facei,                      // label of face
                        mesh.faceNeighbour()[facei],// owner
                        true,                       // face flip
                        newMasterPatches[i],        // patch for face
                        fZone.index(),              // zone for face
                        false,                      // face flip in zone
                        modifiedFace                // modify or add status
                    );
                }

                nModified++;
            }
        }


        // Pass 2. Do other side of zone
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            label zoneFacei = fZone.whichFace(facei);

            if (zoneFacei != -1)
            {
                if (!fZone.flipMap()[zoneFacei])
                {
                    // Use neighbour side of face
                    modifyOrAddFace
                    (
                        meshMod,
                        mesh.faces()[facei].reverseFace(),  // modified face
                        facei,                          // label of face
                        mesh.faceNeighbour()[facei],    // owner
                        true,                           // face flip
                        newSlavePatches[i],             // patch for face
                        fZone.index(),                  // zone for face
                        true,                           // face flip in zone
                        modifiedFace                    // modify or add
                    );
                }
                else
                {
                    // Use owner side of face
                    modifyOrAddFace
                    (
                        meshMod,
                        mesh.faces()[facei],    // modified face
                        facei,                  // label of face
                        mesh.faceOwner()[facei],// owner
                        false,                  // face flip
                        newSlavePatches[i],     // patch for face
                        fZone.index(),          // zone for face
                        true,                   // face flip in zone
                        modifiedFace            // modify or add status
                    );
                }
            }
        }


        // Modify any boundary faces
        // ~~~~~~~~~~~~~~~~~~~~~~~~~

        // Normal boundary:
        // - move to new patch. Might already be back-to-back baffle
        // you want to add cyclic to. Do warn though.
        //
        // Processor boundary:
        // - do not move to cyclic
        // - add normal patches though.

        // For warning once per patch.
        labelHashSet patchWarned;

        forAll(pbm, patchi)
        {
            const polyPatch& pp = pbm[patchi];

            const label newMasterPatchi = newMasterPatches[i];
            const label newSlavePatchi = newSlavePatches[i];

            if
            (
                pp.coupled()
             && (
                    pbm[newMasterPatchi].coupled()
                 || pbm[newSlavePatchi].coupled()
                )
            )
            {
                // Do not allow coupled faces to be moved to different
                // coupled patches.
            }
            else if (pp.coupled() || !internalFacesOnly)
            {
                forAll(pp, i)
                {
                    label facei = pp.start()+i;

                    label zoneFacei = fZone.whichFace(facei);

                    if (zoneFacei != -1)
                    {
//                         if (patchWarned.insert(patchi))
//                         {
//                             WarningInFunction
//                                 << "Found boundary face (in patch "
//                                 << pp.name()
//                                 << ") in faceZone " << fZone.name()
//                                 << " to convert to baffle patches "
//                                 << pbm[newMasterPatchi].name() << "/"
//                                 << pbm[newSlavePatchi].name()
//                                 << endl
//                                 << "    Set internalFacesOnly to true in the"
//                                 << " createBaffles control dictionary if you"
//                                 << " don't wish to convert boundary faces."
//                                 << endl;
//                         }

                        modifyOrAddFace
                        (
                            meshMod,
                            mesh.faces()[facei],        // modified face
                            facei,                      // label of face
                            mesh.faceOwner()[facei],    // owner
                            false,                      // face flip
                            fZone.flipMap()[zoneFacei]
                          ? newSlavePatchi
                          : newMasterPatchi,            // patch for face
                            fZone.index(),              // zone for face
                            fZone.flipMap()[zoneFacei], // face flip in zone
                            modifiedFace                // modify or add
                        );

                        nModified++;
                    }
                }
            }
        }
    }
    return nModified;
}


Foam::label Foam::meshTools::createPatchFaces
(
    const bool internalFacesOnly,
    const polyMesh& mesh,
    const faceZone& fZone,
    const labelList& newPatches,
    polyTopoChange& meshMod,
    PackedBoolList& modifiedFace
)
{
    label nModified = 0;

    forAll(newPatches, i)
    {
        // Pass 1. Do selected side of zone
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            label zoneFacei = fZone.whichFace(facei);

            if (zoneFacei != -1)
            {
                if (!fZone.flipMap()[zoneFacei])
                {
                    // Use owner side of face
                    modifyOrAddFace
                    (
                        meshMod,
                        mesh.faces()[facei],    // modified face
                        facei,                  // label of face
                        mesh.faceOwner()[facei],// owner
                        false,                  // face flip
                        newPatches[i],          // patch for face
                        fZone.index(),          // zone for face
                        false,                  // face flip in zone
                        modifiedFace            // modify or add status
                    );
                }
                else
                {
                    // Use neighbour side of face.
                    // To keep faceZone pointing out of original neighbour
                    // we don't need to set faceFlip since that cell
                    // now becomes the owner
                    modifyOrAddFace
                    (
                        meshMod,
                        mesh.faces()[facei].reverseFace(),  // modified face
                        facei,                      // label of face
                        mesh.faceNeighbour()[facei],// owner
                        true,                       // face flip
                        newPatches[i],              // patch for face
                        fZone.index(),              // zone for face
                        false,                      // face flip in zone
                        modifiedFace                // modify or add status
                    );
                }

                nModified++;
            }
        }

        for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
        {
            label zoneFacei = fZone.whichFace(facei);

            if (zoneFacei != -1)
            {
                // Use owner side of face
                modifyOrAddFace
                (
                    meshMod,
                    mesh.faces()[facei],    // modified face
                    facei,                  // label of face
                    mesh.faceOwner()[facei],// owner
                    false,                  // face flip
                    newPatches[i],          // patch for face
                    fZone.index(),          // zone for face
                    false,                  // face flip in zone
                    modifiedFace            // modify or add status
                );
                nModified++;
            }
        }
    }
    return nModified;
}


Foam::label Foam::meshTools::setRemoveCells
(
    const polyMesh& mesh,
    const labelHashSet& selectedCells,
    const word& patchName,
    polyTopoChange& meshMod,
    const bool keepCells
)
{
    label patchi = mesh.boundaryMesh().findPatchID(patchName);
    labelHashSet cellsToRemove(selectedCells);
    if (keepCells)
    {
        cellsToRemove.clear();
        forAll(mesh.cells(), celli)
        {
            if (!selectedCells.found(celli))
            {
                cellsToRemove.insert(celli);
            }
        }
    }

    removeCells cellRemover(mesh, true);
    labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove.toc()));
    labelList exposedPatchIDs(exposedFaces.size(), patchi);
    cellRemover.setRefinement
    (
        cellsToRemove.toc(),
        exposedFaces,
        exposedPatchIDs,
        meshMod
    );
    return cellsToRemove.size();
}


//- Read and add all fields to the database
void Foam::meshTools::readAndStoreFields(const fvMesh& mesh)
{
    // Get all fields present at the current time
    IOobjectList objects(mesh, mesh.time().timeName());

    readGeoFields<volScalarField>(mesh, objects);
    readGeoFields<volVectorField>(mesh, objects);
    readGeoFields<volSymmTensorField>(mesh, objects);
    readGeoFields<volSphericalTensorField>(mesh, objects);
    readGeoFields<volTensorField>(mesh, objects);

    readGeoFields<volScalarField::Internal>(mesh, objects);
    readGeoFields<volVectorField::Internal>(mesh, objects);
    readGeoFields<volSymmTensorField::Internal>(mesh, objects);
    readGeoFields<volSphericalTensorField::Internal>(mesh, objects);
    readGeoFields<volTensorField::Internal>(mesh, objects);

    readGeoFields<surfaceScalarField>(mesh, objects);
    readGeoFields<surfaceVectorField>(mesh, objects);
    readGeoFields<surfaceSymmTensorField>(mesh, objects);
    readGeoFields<surfaceSphericalTensorField>(mesh, objects);
    readGeoFields<surfaceTensorField>(mesh, objects);

    readPointFields<pointScalarField>(mesh, objects);
    readPointFields<pointVectorField>(mesh, objects);
    readPointFields<pointSymmTensorField>(mesh, objects);
    readPointFields<pointSphericalTensorField>(mesh, objects);
    readPointFields<pointTensorField>(mesh, objects);
}


void Foam::meshTools::filterPatches
(
    polyMesh& mesh,
    const HashSet<word>& addedPatchNames
)
{
    // Remove any zero-sized ones. Assumes
    // - processor patches are already only there if needed
    // - all other patches are available on all processors
    // - but coupled ones might still be needed, even if zero-size
    //   (e.g. processorCyclic)
    // See also logic in createPatch.
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    labelList oldToNew(pbm.size(), -1);
    label newPatchi = 0;
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (!isA<processorPolyPatch>(pp))
        {
            if
            (
                isA<coupledPolyPatch>(pp)
             || returnReduce(pp.size(), sumOp<label>())
             || addedPatchNames.found(pp.name())
            )
            {
                // Coupled (and unknown size) or uncoupled and used
                oldToNew[patchi] = newPatchi++;
            }
        }
    }

    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            oldToNew[patchi] = newPatchi++;
        }
    }


    const label nKeepPatches = newPatchi;

    // Shuffle unused ones to end
    if (nKeepPatches != pbm.size())
    {
        Info<< endl
            << "Removing zero-sized patches:" << endl << incrIndent;

        forAll(oldToNew, patchi)
        {
            if (oldToNew[patchi] == -1)
            {
                Info<< indent << pbm[patchi].name()
                    << " type " << pbm[patchi].type()
                    << " at position " << patchi << endl;
                oldToNew[patchi] = newPatchi++;
            }
        }
        Info<< decrIndent;

        mesh.reorderPatches(oldToNew, true);
        Info<< endl;
    }
}

// ************************************************************************* //
