/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meshTools::setFaceInfo
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
//                 WarningInFunction
//                     << pp.name() << " and " << pbm[newMasterPatchi].name()
//                     << "/" << pbm[newSlavePatchi].name()
//                     << " are both coupled "
//                     << " so these shared faces will not be modified." << nl
//                     << "If this is a processor patch, try changing the "
//                     << "distribution method." << endl;
            }
            else if (pp.coupled() || !internalFacesOnly)
            {
                forAll(pp, i)
                {
                    label facei = pp.start()+i;

                    label zoneFacei = fZone.whichFace(facei);

                    if (zoneFacei != -1)
                    {
                        if (patchWarned.insert(patchi))
                        {
                            WarningInFunction
                                << "Found boundary face (in patch "
                                << pp.name()
                                << ") in faceZone " << fZone.name()
                                << " to convert to baffle patches "
                                << pbm[newMasterPatchi].name() << "/"
                                << pbm[newSlavePatchi].name()
                                << endl
                                << "    Set internalFacesOnly to true in the"
                                << " createBaffles control dictionary if you"
                                << " don't wish to convert boundary faces."
                                << endl;
                        }

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
    fvMesh& mesh,
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

        fvMeshTools::reorderPatches(mesh, oldToNew, nKeepPatches, true);
        Info<< endl;
    }
}

// ************************************************************************* //
