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

#include "cellRemovalFvMesh.H"
#include "regionSplit.H"
#include "mapPolyMesh.H"
#include "removeCells.H"
#include "polyTopoChange.H"
#include "volMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellRemovalFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, cellRemovalFvMesh, IOobject);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::cellRemovalFvMesh::cellRemovalFvMesh
(
    const IOobject& io
)
:
    dynamicFvMesh(io),
    dict_(dynamicMeshDict().optionalSubDict(type() + "Coeffs")),
    removeDeadCells_
    (
        dict_.lookupOrDefault<Switch>
        (
            "removeDeadCells",
            true
        )
    ),
    lawPtr_
    (
        cellRemovalLaw::New
        (
            "law",
            *this,
            dict_.subDict("law")
        )
     )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellRemovalFvMesh::~cellRemovalFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cellRemovalFvMesh::update()
{
    // Check if there are cells to remove
    const labelField cellsToRemove(lawPtr_->cellsToRemove());

    const label nCellsToRemove =
        returnReduce(cellsToRemove.size(), sumOp<int>());

    if (nCellsToRemove)
    {
        // Exposed faces will be inserted into the open patch

        Info<< nl << "There are " << nCellsToRemove
            << " cells to be removed" << endl;

        // Find faces that will be exposed
        // These faces will become boundary faces

        removeCells cellRemover(*this);
        const labelList facesToExpose
        (
            cellRemover.getExposedFaces(cellsToRemove)
        );

        Info<< "There are " << facesToExpose.size()
            << " internal faces that will be exposed" << endl;

        // Clear the directTopoChanger before giving it to the cell remover
        clearOut();

        // Set actions in cell remover
        polyTopoChange meshMod(*this);
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
        Info<< "Performing mesh change" << endl;
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, true);

        // Update cell remover map
        cellRemover.updateMesh(map());

        // Update mesh fields e.g. U, sigma, etc.
        updateMesh(map);

        // Update fields on newly exposed faces
        Info<< "Updating field values on newly exposed faces" << endl;
        updateVolFieldsExposedFaces<scalar>(map, facesToExpose);
        updateVolFieldsExposedFaces<vector>(map, facesToExpose);
        updateVolFieldsExposedFaces<tensor>(map, facesToExpose);
        updateVolFieldsExposedFaces<symmTensor>(map, facesToExpose);
//         updateVolFieldsExposedFaces<diagTensor>(map, facesToExpose);
        updateVolFieldsExposedFaces<sphericalTensor>(map, facesToExpose);

        // Move mesh (since morphing does not do this)
        if (map().hasMotionPoints())
        {
            fvMesh::movePoints(map().preMotionPoints());
        }

        // Remove any dead cells
        if (removeDeadCells_)
        {
            Info<< "Checking for dead cells" << endl;

            regionSplit regions(*this);

            // Make nCellsInRegion list

            labelList nCellsInRegion(regions.nRegions(), 0);

            forAll(regions, cellI)
            {
                nCellsInRegion[regions[cellI]]++;
            }

            // Find regions with less cells than the dead-cell limit
            Info<< "The number of ncells is " << nCells()<< endl;

            const label deadCellLimit = 0.1*nCells();
            labelHashSet deadCellsSet;

            forAll(regions, cellI)
            {
                if (nCellsInRegion[regions[cellI]] < deadCellLimit)
                {
                    deadCellsSet.insert(cellI);
                }
            }

            const labelList deadCells = deadCellsSet.toc();

            label ndeadCells = deadCells.size();
            reduce(ndeadCells, maxOp<label>());

            if (ndeadCells)
            {
                // Remove dead cells from the mesh

                // Clear the directTopoChanger before giving it to the cell
                // remover
                clearOut();

                // Set actions in cell remover
                cellRemover.setRefinement
                (
                    deadCells,
                    labelList(0),    // no faces to expose
                    labelList(0),    // open patch not needed
                    meshMod
                );

                // Change the mesh
                Info<< "    Removing " << deadCells.size() << " dead cells"
                    << endl;

                autoPtr<mapPolyMesh> dmap =
                    meshMod.changeMesh(*this, true);

                // Update cell remover map
                cellRemover.updateMesh(dmap());

                // Update mesh fields e.g. U, sigma, etc.
                updateMesh(dmap);

                // Move mesh (since morphing does not do this)
                if (dmap().hasMotionPoints())
                {
                    fvMesh::movePoints(dmap().preMotionPoints());
                }
            }
            else
            {
                Info<< "    There are no dead cells" << endl;
            }
        }
    }

    if (time().outputTime())
    {
        write();
    }

    return bool(nCellsToRemove > 0);
}


// ************************************************************************* //
