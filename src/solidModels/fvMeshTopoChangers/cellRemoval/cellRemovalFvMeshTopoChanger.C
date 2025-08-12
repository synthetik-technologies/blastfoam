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

#include "cellRemovalFvMeshTopoChanger.H"
#include "regionSplit.H"
#include "polyTopoChangeMap.H"
#include "removeCells.H"
#include "polyTopoChange.H"
#include "volMesh.H"
#include "lookupSolidModel.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{
    defineTypeNameAndDebug(cellRemoval, 0);
    addToRunTimeSelectionTable(fvMeshTopoChanger, cellRemoval, fvMesh);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::cellRemoval::cellRemoval
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    fvMeshTopoChanger(mesh),
    removeDeadCells_
    (
        dict.lookupOrDefault<Switch>
        (
            "removeDeadCells",
            false
        )
    ),
    lawPtr_(cellRemovalLaw::New("law", mesh, dict)),
    saveSubset_(false),
    subsetter_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::cellRemoval::~cellRemoval()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshTopoChangers::cellRemoval::update()
{
    fvMesh& mesh = this->mesh();

    // Check if there are cells to remove
    labelList cellsToRemove(lawPtr_->cellsToRemove().toc());

    const label nCellsToRemove =
        returnReduce(cellsToRemove.size(), sumOp<label>());

    subsetter_.clear();

    if (!nCellsToRemove)
    {
        return false;
    }

    Info<< "Marked " << nCellsToRemove << " cells for removal" << endl;
    if (removeDeadCells_)
    {
        const label nDeadCells = addDeadCells(cellsToRemove);
        Info<< "Marked " << returnReduce(nDeadCells, sumOp<label>())
            << " additional dead cells for removal" << endl;
    }

    // Save the removed region to a subsetMesh
    if (saveSubset_)
    {
        // Create the subsetter
        subsetter_.set(new fvMeshSubset(mesh));

        // Subset the mesh
        subsetter_->setLargeCellSubset(cellsToRemove);

        // Rename the mesh
        subsetter_->subMesh().polyMesh::rename(mesh.name() + "_removed");

        // Lookup the solid model
        const solidModel& solid = lookupSolidModel(mesh);

        // If the mesh is not moving, move the subset mesh to the deformed
        // geometry
        if (!solid.movingMesh())
        {
            subsetter_->subMesh().movePoints
            (
                subsetter_->subMesh().points()
              + pointField(solid.pointD(), subsetter_->pointMap())
            );
        }

        #define saveGeoFieldTypes(Type, Patch, Mesh) \
            saveGeoFields<Type, Patch, Mesh>(mesh);

        FOR_ALL_FIELD_TYPES(saveGeoFieldTypes, fvPatchField, volMesh)
        #undef saveGeoFieldTypes
    }

    const label nOldCells = returnReduce(mesh.nCells(), sumOp<label>());

    // Exposed faces will be inserted into the open patch
    Info<< nl << "Selected " << nCellsToRemove
        << " cells to remove" << endl;

    // Find faces that will be exposed
    // These faces will become boundary faces

    removeCells cellRemover(mesh);
    const labelList facesToExpose
    (
        cellRemover.getExposedFaces(cellsToRemove)
    );

    DebugInfo
        << "There are " << returnReduce(facesToExpose.size(), sumOp<label>())
        << " internal faces that will be exposed" << endl;

    // Set actions in cell remover
    polyTopoChange meshMod(mesh);
    cellRemover.setRefinement
    (
        cellsToRemove,
        facesToExpose,
        labelList
        (
            facesToExpose.size(),
            lawPtr_->exposedFacesPatchID()
        ),
        meshMod
    );

    // Change the mesh
    DebugInfo<< "Performing mesh change" << endl;
    autoPtr<polyTopoChangeMap> map = meshMod.changeMesh(mesh);

    // Update mesh fields e.g. U, sigma, etc.
    mesh.topoChange(map);

    labelList pfMap(patchFaceMap(map, facesToExpose));

    // Update fields on newly exposed faces
    DebugInfo<< "Updating field values on newly exposed faces" << endl;
    updateVolFieldsExposedFaces<scalar>(map, facesToExpose, pfMap);
    updateVolFieldsExposedFaces<vector>(map, facesToExpose, pfMap);
    updateVolFieldsExposedFaces<tensor>(map, facesToExpose, pfMap);
    updateVolFieldsExposedFaces<symmTensor>(map, facesToExpose, pfMap);
    updateVolFieldsExposedFaces<sphericalTensor>(map, facesToExpose, pfMap);

    Info<< "Changed from " << nOldCells << " cells to "
            << returnReduce(mesh.nCells(), sumOp<label>()) << " cells"<<endl;

    return nCellsToRemove > 0;
}


Foam::label Foam::fvMeshTopoChangers::cellRemoval::addDeadCells
(
    labelList& cellsToRemove
)
{
    const label nOldCellsToRemove = cellsToRemove.size();

    // Lookup the solid model
    const fvMesh& mesh = this->mesh();
    const solidModel& solid = lookupSolidModel(mesh);
    const volVectorField& D = solid.solutionD();
    const cellList& cells = mesh.cells();
    const labelList& owner = mesh.faceOwner();
    const labelList& neighbour = mesh.faceNeighbour();
    PackedBoolList markedCells(mesh.nCells(), true);
    PackedBoolList markedFaces(mesh.nFaces(), true);
    PackedBoolList visitedCells(mesh.nCells());
    PackedBoolList visitedFaces(mesh.nFaces());

    // Faces to visit on the next pass
    labelHashSet newFacesToVisit;

    // Faces to visit this pass
    labelHashSet facesToVisit;

    //- Mark boundary faces that fix values
    forAll(D.boundaryField(), patchi)
    {
        if (D.boundaryField()[patchi].fixesValue())
        {
            const polyPatch& patch = mesh.boundaryMesh()[patchi];
            const labelList& faceCells = patch.faceCells();
            forAll(patch, fi)
            {
                const label celli = faceCells[fi];
                const cell& c = cells[celli];
                forAll(c, fj)
                {
                    const label facej = c[fj];
                    markedFaces.unset(facej);
                    visitedFaces.set(facej);
                    facesToVisit.set(facej);
                }
                markedCells.unset(celli);
                visitedCells.set(celli);
            }
        }
    }

    //- Mark cell faces with fixed values
    forAll(solid.setCellDisps().cellIDs(), ci)
    {
        const label celli = solid.setCellDisps().cellIDs()[ci];
        if (!visitedFaces.get(celli))
        {
            const cell& c = cells[celli];
            forAll(c, fi)
            {
                const label facei = c[fi];
                markedFaces.unset(facei);
                visitedFaces.set(facei);
                facesToVisit.set(facei);
            }
            markedCells.unset(celli);
            visitedCells.set(celli);
        }
    }

    //- Unmark faces with cells to be removed
    forAll(cellsToRemove, ci)
    {
        const label celli = cellsToRemove[ci];
        const cell& c = cells[celli];
        forAll(c, fi)
        {
            const label facei = c[fi];
            markedFaces.set(facei);
            visitedFaces.set(facei);
            facesToVisit.erase(facei);
        }
        markedCells.set(celli);
        visitedCells.set(celli);
    }

    //- Synchronize boundary faces across coupled patches
    syncVisitedFaces(visitedFaces, markedFaces, newFacesToVisit);

    // Loop until no more valid faces to visit
    do
    {
        // Clear the list of faces to visit
        newFacesToVisit.clear();

        // Only loop thorough faces that were added in the previous iteration
        forAllConstIter(labelHashSet, facesToVisit, iter)
        {
            const label facei = iter.key();

            // Check if there is a owner/neighbour cell that has not yet been
            // visited
            label celli = -1;
            if (!visitedCells.get(owner[facei]))
            {
                celli = owner[facei];
            }
            else if
            (
                facei < mesh.nInternalFaces()
             && !visitedCells.get(neighbour[facei])
            )
            {
                celli = neighbour[facei];
            }

            // The owner/neighbour has not been added yet so visit and mark
            // all of its faces if not visited yet
            if (celli >= 0)
            {
                const cell& c = cells[celli];
                forAll(c, fj)
                {
                    const label facej = c[fj];

                    // Only unmark if not yet visited
                    if (visitedFaces.set(facej))
                    {
                        newFacesToVisit.set(facej);
                        markedFaces.unset(facej);
                    }
                }

                // Mark the cell as visited and unmark it for removal since we
                // were able to walk here
                visitedCells.set(celli);
                markedCells.unset(celli);
            }
        }

        //- Synchronize boundary faces across coupled patches
        syncVisitedFaces(visitedFaces, markedFaces, newFacesToVisit);

        // Done with this iteration so transfer the newly visited faces to the
        // list of faces to visit next iteration
        facesToVisit.transfer(newFacesToVisit);

    } while (returnReduce(facesToVisit.size(), sumOp<label>()));

    cellsToRemove.resize(markedCells.count());
    label ci = 0;
    forAll(markedCells, celli)
    {
        if (markedCells.get(celli))
        {
            cellsToRemove[ci++] = celli;
        }
    }
    return markedCells.count() - nOldCellsToRemove;
}


void Foam::fvMeshTopoChangers::cellRemoval::syncVisitedFaces
(
    PackedBoolList& visited,
    PackedBoolList& marked,
    labelHashSet& facesToVisit
) const
{
    if (Pstream::parRun())
    {
        const label nInternalFaces = mesh().nInternalFaces();
        const polyBoundaryMesh& bmesh = mesh().boundaryMesh();

        // Current state of boundary faces
        // 0: Not visited
        // 1: Visited and marked
        // 2: Visited and unmarked
        labelList boundaryState(mesh().nFaces() - nInternalFaces, 0);
        forAll(bmesh, patchi)
        {
            const polyPatch& patch = bmesh[patchi];
            const label start = patch.start();
            forAll(patch, fi)
            {
                const label facei = start + fi;
                if (visited.get(facei))
                {
                    boundaryState[facei - nInternalFaces] =
                        marked.get(facei) ? 1 : 2;
                }
            }
        }

        syncTools::syncBoundaryFaceList
        (
            mesh(),
            boundaryState,
            maxEqOp<label>()
        );

        forAll(bmesh, patchi)
        {
            const polyPatch& patch = bmesh[patchi];
            const label start = patch.start();
            forAll(patch, fi)
            {
                const label facei = start + fi;
                const label bfacei = facei - nInternalFaces;
                if (boundaryState[bfacei] > 0)
                {
                    if (visited.set(facei))
                    {
                        facesToVisit.insert(facei);
                        if (boundaryState[bfacei] > 1)
                        {
                            marked.unset(facei);
                        }
                    }
                }
            }
        }
    }
}


Foam::labelList Foam::fvMeshTopoChangers::cellRemoval::patchFaceMap
(
    const polyTopoChangeMap& mpm,
    const labelList& exposedFaces
) const
{
    labelList map(exposedFaces.size(), -1);

    // Get reverse face map
    const labelList& revFaceMap = mpm.reverseFaceMap();

    // Set the patch index of newly exposed faces using the old index
    forAll(exposedFaces, fi)
    {
        // Get new face ID
        label newFaceID = revFaceMap[exposedFaces[fi]];

        // Find the patch ID
        const label patchID = mesh().boundaryMesh().whichPatch(newFaceID);

        if (patchID == -1)
        {
            FatalErrorInFunction
                << "exposed face is not on the boundary!? What's going on?"
                << abort(FatalError);
        }

        map[fi] = patchID;
    }

    return map;
}


void Foam::fvMeshTopoChangers::cellRemoval::topoChange
(
    const polyTopoChangeMap& map
)
{}


void Foam::fvMeshTopoChangers::cellRemoval::mapMesh(const polyMeshMap& map)
{}


void Foam::fvMeshTopoChangers::cellRemoval::distribute
(
    const polyDistributionMap& map
)
{}

// ************************************************************************* //
