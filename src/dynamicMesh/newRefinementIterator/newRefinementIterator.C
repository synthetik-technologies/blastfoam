/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "newRefinementIterator.H"
#include "polyMesh.H"
#include "Time.H"
#include "refineCell.H"
#include "undoableMeshCutter.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "newCellCuts.H"
#include "OFstream.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(newRefinementIterator, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::newRefinementIterator::newRefinementIterator
(
    polyMesh& mesh,
    undoableMeshCutter& meshRefiner,
    const cellLooper& cellWalker
)
:
    edgeVertex(mesh),
    mesh_(mesh),
    meshRefiner_(meshRefiner),
    cellWalker_(cellWalker)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::newRefinementIterator::~newRefinementIterator()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Map<Foam::label> Foam::newRefinementIterator::setRefinement
(
    const List<refineCell>& refCells
)
{
    Map<label> addedCells(2*refCells.size());

    Time& runTime = const_cast<Time&>(mesh_.time());

    label nRefCells = refCells.size();

    label oldRefCells = -1;

    // Operate on copy.
    List<refineCell> currentRefCells(refCells);

    bool stop = false;

    do
    {
        polyTopoChange meshMod(mesh_);

        if (debug)
        {
            Pout<< "newRefinementIterator : refining "
                << currentRefCells.size() << " cells." << endl;
        }

        // Determine cut pattern.
        cellCuts cuts(mesh_, cellWalker_, currentRefCells);

        label nCuts = cuts.nLoops();
        reduce(nCuts, sumOp<label>());

        if (nCuts == 0)
        {
            if (debug)
            {
                Pout<< "newRefinementIterator : exiting iteration since no valid"
                    << " loops found for " << currentRefCells.size()
                    << " cells" << endl;


                fileName cutsFile("failedCuts_" + runTime.timeName() + ".obj");

                Pout<< "Writing cuts for time " <<  runTime.timeName()
                    << " to " << cutsFile << endl;

                OFstream cutsStream(cutsFile);


                labelList refCellsDebug(currentRefCells.size());
                forAll(currentRefCells, i)
                {
                    refCellsDebug[i] = currentRefCells[i].cellNo();
                }
                meshTools::writeOBJ
                (
                    cutsStream,
                    mesh().cells(),
                    mesh().faces(),
                    mesh().points(),
                    refCellsDebug
                );
            }

            break;
        }

        if (debug)
        {
            fileName cutsFile("cuts_" + runTime.timeName() + ".obj");

            Pout<< "Writing cuts for time " <<  runTime.timeName()
                << " to " << cutsFile << endl;

            OFstream cutsStream(cutsFile);
            cuts.writeOBJ(cutsStream);
        }


        // Insert mesh refinement into polyTopoChange.
        meshRefiner_.setRefinement(cuts, meshMod);


        //
        // Do all changes
        //

        autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh
        (
            mesh_,
            false
        );

        // Move mesh (since morphing does not do this)
        if (morphMap().hasMotionPoints())
        {
            mesh_.movePoints(morphMap().preMotionPoints());
        }

        // Update stored refinement pattern
        mesh_.updateMesh(morphMap());
        meshRefiner_.updateMesh(morphMap());

        // Update currentRefCells for new cell numbers. Use helper function
        // in meshCutter class.
        updateLabels
        (
            morphMap->reverseCellMap(),
            currentRefCells
        );

        // Update addedCells for new cell numbers
        updateLabels
        (
            morphMap->reverseCellMap(),
            addedCells
        );

        // Get all added cells from cellCutter (already in new numbering
        // from meshRefiner.updateMesh call) and add to global list of added
        const Map<label>& addedNow = meshRefiner_.addedCells();

        forAllConstIter(Map<label>, addedNow, iter)
        {
            if (!addedCells.insert(iter.key(), iter()))
            {
                FatalErrorInFunction
                    << "Master cell " << iter.key()
                    << " already has been refined" << endl
                    << "Added cell:" << iter() << abort(FatalError);
            }
        }


        // Get failed refinement in new cell numbering and reconstruct input
        // to the meshRefiner. Is done by removing all refined cells from
        // current list of cells to refine.

        // Update refCells for new cell numbers.
        updateLabels
        (
            morphMap->reverseCellMap(),
            currentRefCells
        );

        // Pack refCells acc. to refined status
        nRefCells = 0;

        forAll(currentRefCells, refI)
        {
            const refineCell& refCell = currentRefCells[refI];

            if (!addedNow.found(refCell.cellNo()))
            {
                if (nRefCells != refI)
                {
                    currentRefCells[nRefCells++] =
                        refineCell
                        (
                            refCell.cellNo(),
                            refCell.direction()
                        );
                }
            }
        }

        oldRefCells = currentRefCells.size();

        currentRefCells.setSize(nRefCells);

        if (debug)
        {
            Pout<< endl;
        }

        // Stop only if all finished or all can't refine any further.
        stop = (nRefCells == 0) || (nRefCells == oldRefCells);
        reduce(stop, andOp<bool>());
    }
    while (!stop);


    if (returnReduce((nRefCells == oldRefCells), andOp<bool>()))
    {
        WarningInFunction
            << "stopped refining."
            << "Did not manage to refine a single cell" << endl
            << "Wanted :" << oldRefCells << endl;
    }

    return addedCells;
}


Foam::Map<Foam::label> Foam::newRefinementIterator::setRefinement
(
    const List<label>& refCells,
    const PtrList<plane>& refPlanes
)
{
    Map<label> addedCells(2*refCells.size());

    Time& runTime = const_cast<Time&>(mesh_.time());

    label nRefCells = refCells.size();

    label oldRefCells = -1;

    // Operate on copy.
    List<label> currentRefCells(refCells);
    PtrList<plane> currentRefPlanes(refPlanes.size());
    forAll(currentRefPlanes, i)
    {
        currentRefPlanes.set(i, new plane(refPlanes[i]));
    }

    bool stop = false;

    do
    {
        polyTopoChange meshMod(mesh_);

        if (debug)
        {
            Pout<< "newRefinementIterator : refining "
                << currentRefCells.size() << " cells." << endl;
        }

        // Determine cut pattern.
        newCellCuts cuts(mesh_, cellWalker_, currentRefCells, currentRefPlanes);

        label nCuts = cuts.nLoops();
        reduce(nCuts, sumOp<label>());

        if (nCuts == 0)
        {
            if (debug)
            {
                Pout<< "newRefinementIterator : exiting iteration since no valid"
                    << " loops found for " << currentRefCells.size()
                    << " cells" << endl;


                fileName cutsFile("failedCuts_" + runTime.timeName() + ".obj");

                Pout<< "Writing cuts for time " <<  runTime.timeName()
                    << " to " << cutsFile << endl;

                OFstream cutsStream(cutsFile);


                labelList refCellsDebug(currentRefCells.size());
                forAll(currentRefCells, i)
                {
                    refCellsDebug[i] = currentRefCells[i];
                }
                meshTools::writeOBJ
                (
                    cutsStream,
                    mesh().cells(),
                    mesh().faces(),
                    mesh().points(),
                    refCellsDebug
                );
            }

            break;
        }

        if (debug)
        {
            fileName cutsFile("cuts_" + runTime.timeName() + ".obj");

            Pout<< "Writing cuts for time " <<  runTime.timeName()
                << " to " << cutsFile << endl;

            OFstream cutsStream(cutsFile);
            cuts.writeOBJ(cutsStream);
        }


        // Insert mesh refinement into polyTopoChange.
        meshRefiner_.setRefinement(cuts, meshMod);


        //
        // Do all changes
        //

        autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh
        (
            mesh_,
            false
        );

        // Move mesh (since morphing does not do this)
        if (morphMap().hasMotionPoints())
        {
            mesh_.movePoints(morphMap().preMotionPoints());
        }

        // Update stored refinement pattern
        mesh_.updateMesh(morphMap());
        meshRefiner_.updateMesh(morphMap());

        // Update addedCells for new cell numbers
        updateLabels
        (
            morphMap->reverseCellMap(),
            addedCells
        );
        {
            const labelList& map = morphMap->reverseCellMap();
            label newRefI = 0;

            forAll(currentRefPlanes, refI)
            {
                const label oldCelli = currentRefCells[refI];
                const plane& refPlane = currentRefPlanes[refI];

                label newCelli = map[oldCelli];

                if (newCelli != -1)
                {
                    currentRefCells[newRefI] = newCelli;
                    currentRefPlanes[newRefI++] = plane(refPlane);
                }
            }
            currentRefCells.setSize(newRefI);
            currentRefPlanes.setSize(newRefI);
        }

        // Get all added cells from cellCutter (already in new numbering
        // from meshRefiner.updateMesh call) and add to global list of added
        const Map<label>& addedNow = meshRefiner_.addedCells();

        forAllConstIter(Map<label>, addedNow, iter)
        {
            if (!addedCells.insert(iter.key(), iter()))
            {
                FatalErrorInFunction
                    << "Master cell " << iter.key()
                    << " already has been refined" << endl
                    << "Added cell:" << iter() << abort(FatalError);
            }
        }


        // Get failed refinement in new cell numbering and reconstruct input
        // to the meshRefiner. Is done by removing all refined cells from
        // current list of cells to refine.

        // Update refCells for new cell numbers.
        {
            const labelList& map = morphMap->reverseCellMap();
            label newRefI = 0;

            forAll(currentRefPlanes, refI)
            {
                const label oldCelli = currentRefCells[refI];
                const plane& refPlane = currentRefPlanes[refI];

                label newCelli = map[oldCelli];

                if (newCelli != -1)
                {
                    currentRefCells[newRefI] = newCelli;
                    currentRefPlanes[newRefI++] = plane(refPlane);
                }
            }
            currentRefCells.setSize(newRefI);
            currentRefPlanes.setSize(newRefI);
        }

        // Pack refCells acc. to refined status
        nRefCells = 0;

        forAll(currentRefCells, refI)
        {
            const label refCell = currentRefCells[refI];
            const plane& refPlane = currentRefPlanes[refI];

            if (!addedNow.found(refCell))
            {
                if (nRefCells != refI)
                {
                    currentRefCells[nRefCells++] = refCell;
                    currentRefPlanes[nRefCells++] = refPlane;
                }
            }
        }

        oldRefCells = currentRefCells.size();

        currentRefCells.setSize(nRefCells);
        currentRefPlanes.setSize(nRefCells);

        if (debug)
        {
            Pout<< endl;
        }

        // Stop only if all finished or all can't refine any further.
        stop = (nRefCells == 0) || (nRefCells == oldRefCells);
        reduce(stop, andOp<bool>());
    }
    while (!stop);


    if (returnReduce((nRefCells == oldRefCells), andOp<bool>()))
    {
        WarningInFunction
            << "stopped refining."
            << "Did not manage to refine a single cell" << endl
            << "Wanted :" << oldRefCells << endl;
    }

    return addedCells;
}


// ************************************************************************* //
