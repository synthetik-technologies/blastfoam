/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "regionsToCell.H"
#include "regionSplit.H"
#include "emptyPolyPatch.H"
#include "cellSet.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionsToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, regionsToCell, word);
    addToRunTimeSelectionTable(topoSetSource, regionsToCell, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetSource,
        regionsToCell,
        word,
        regions
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetSource,
        regionsToCell,
        istream,
        regions
    );
}


Foam::topoSetSource::addToUsageTable Foam::regionsToCell::usage_
(
    regionsToCell::typeName,
    "\n    Usage: regionsToCell subCellSet (pt0 .. ptn) nErode\n\n"
    "    Select all cells in the connected region containing"
    " points (pt0..ptn).\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::regionsToCell::markRegionFaces
(
    const boolList& selectedCell,
    boolList& regionFace
) const
{
    // Internal faces
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();
    forAll(faceNeighbour, faceI)
    {
        if
        (
            selectedCell[faceOwner[faceI]]
         != selectedCell[faceNeighbour[faceI]]
        )
        {
            regionFace[faceI] = true;
        }
    }

    // Swap neighbour selectedCell state
    boolList nbrSelected;
    syncTools::swapBoundaryCellList(mesh_, selectedCell, nbrSelected);

    // Boundary faces
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];
        const labelUList& faceCells = pp.faceCells();
        forAll(faceCells, i)
        {
            label faceI = pp.start()+i;
            label bFaceI = faceI-mesh_.nInternalFaces();
            if
            (
                selectedCell[faceCells[i]]
             != selectedCell[nbrSelected[bFaceI]]
            )
            {
                regionFace[faceI] = true;
            }
        }
    }
}


Foam::boolList Foam::regionsToCell::findRegions
(
    const bool verbose,
    const boolList& selectedCell,
    const regionSplit& cellRegion
) const
{
    boolList keepRegion(cellRegion.nRegions(), false);

    for (const point& insidePt : insidePoints_)
    {
        // Find the region containing the insidePoint

        //label cellI = mesh_.findCell(insidePt);
        label cellI = -1;
        forAll(selectedCell, index)
        {
            if
            (
                selectedCell[index]
             && mesh_.pointInCell(insidePt, index, polyMesh::CELL_TETS)
            )
            {
                cellI = index;
                break;
            }
        }

        label keepRegionI = -1;
        label keepProcI = -1;
        if (cellI != -1)
        {
            keepRegionI = cellRegion[cellI];
            keepProcI = Pstream::myProcNo();
        }
        reduce(keepRegionI, maxOp<label>());
        keepRegion[keepRegionI] = true;

        reduce(keepProcI, maxOp<label>());

        if (keepProcI == -1)
        {
            FatalErrorInFunction
                << "Did not find " << insidePt
                << " in mesh." << " Mesh bounds are " << mesh_.bounds()
                << exit(FatalError);
        }

        if (verbose)
        {
            Info<< "    Found location " << insidePt
                << " in cell " << cellI << " on processor " << keepProcI
                << " in global region " << keepRegionI
                << " out of " << cellRegion.nRegions() << " regions." << endl;
        }
    }

    return keepRegion;
}


void Foam::regionsToCell::unselectOutsideRegions
(
    boolList& selectedCell
) const
{
    // Determine faces on the edge of selectedCell
    boolList blockedFace(mesh_.nFaces(), false);
    markRegionFaces(selectedCell, blockedFace);

    // Determine regions
    regionSplit cellRegion(mesh_, blockedFace);

    // Determine regions containing insidePoints_
    boolList keepRegion(findRegions(debug, selectedCell, cellRegion));

    // Go back to bool per cell
    forAll(cellRegion, cellI)
    {
        if (!keepRegion[cellRegion[cellI]])
        {
            selectedCell[cellI] = false;
        }
    }
}


void Foam::regionsToCell::shrinkRegions
(
    boolList& selectedCell
) const
{
    // Select points on unselected cells and boundary
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    boolList boundaryPoint(mesh_.nPoints(), false);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    for (const polyPatch& pp : pbm)
    {
        if (!pp.coupled() && !isA<emptyPolyPatch>(pp))
        {
            forAll(pp, i)
            {
                const face& f = pp[i];
                forAll(f, fp)
                {
                    boundaryPoint[f[fp]] = true;
                }
            }
        }
    }

    forAll(selectedCell, cellI)
    {
        if (!selectedCell[cellI])
        {
            const labelList& cPoints = mesh_.cellPoints(cellI);
            forAll(cPoints, i)
            {
                boundaryPoint[cPoints[i]] = true;
            }
        }
    }

    syncTools::syncPointList(mesh_, boundaryPoint, orEqOp<bool>(), false);


    // Select all cells using these points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    label nChanged = 0;
    forAll(boundaryPoint, pointI)
    {
        if (boundaryPoint[pointI])
        {
            const labelList& pCells = mesh_.pointCells(pointI);
            forAll(pCells, i)
            {
                label cellI = pCells[i];
                if (selectedCell[cellI])
                {
                    selectedCell[cellI] = false;
                    nChanged++;
                }
            }
        }
    }
    Info<< "    Eroded " << returnReduce(nChanged, sumOp<label>())
        << " cells." << endl;
}


void Foam::regionsToCell::erode
(
    boolList& selectedCell
) const
{
    //Info<< "Entering shrinkRegions:" << count(selectedCell) << endl;
    //generateField("selectedCell_before", selectedCell)().write();

    // Now erode and see which regions get disconnected
    boolList shrunkSelectedCell(selectedCell);

    for (label iter = 0; iter < nErode_; iter++)
    {
        shrinkRegions(shrunkSelectedCell);
    }

    //Info<< "After shrinking:" << count(shrunkSelectedCell) << endl;
    //generateField("shrunkSelectedCell", shrunkSelectedCell)().write();



    // Determine faces on the edge of shrunkSelectedCell
    boolList blockedFace(mesh_.nFaces(), false);
    markRegionFaces(shrunkSelectedCell, blockedFace);

    // Find disconnected regions
    regionSplit cellRegion(mesh_, blockedFace);

    // Determine regions containing insidePoints
    boolList keepRegion(findRegions(debug, shrunkSelectedCell, cellRegion));


    // Extract cells in regions that are not to be kept.
    boolList removeCell(mesh_.nCells(), false);
    forAll(cellRegion, cellI)
    {
        if (shrunkSelectedCell[cellI] && !keepRegion[cellRegion[cellI]])
        {
            removeCell[cellI] = true;
        }
    }

    //Info<< "removeCell before:" << count(removeCell) << endl;
    //generateField("removeCell_before", removeCell)().write();



    // Grow removeCell
    for (label iter = 0; iter < nErode_; iter++)
    {
        // Grow selected cell in regions that are not for keeping
        boolList boundaryPoint(mesh_.nPoints(), false);
        forAll(removeCell, cellI)
        {
            if (removeCell[cellI])
            {
                const labelList& cPoints = mesh_.cellPoints(cellI);
                forAll(cPoints, i)
                {
                    boundaryPoint[cPoints[i]] = true;
                }
            }
        }
        syncTools::syncPointList(mesh_, boundaryPoint, orEqOp<bool>(), false);

        // Select all cells using these points

        label nChanged = 0;
        forAll(boundaryPoint, pointI)
        {
            if (boundaryPoint[pointI])
            {
                const labelList& pCells = mesh_.pointCells(pointI);
                forAll(pCells, i)
                {
                    label cellI = pCells[i];
                    if (!removeCell[cellI])
                    {
                        removeCell[cellI] = true;
                        nChanged++;
                    }
                }
            }
        }
    }

    //Info<< "removeCell after:" << count(removeCell) << endl;
    //generateField("removeCell_after", removeCell)().write();


    // Unmark removeCell
    forAll(removeCell, cellI)
    {
        if (removeCell[cellI])
        {
            selectedCell[cellI] = false;
        }
    }
}


void Foam::regionsToCell::combine(topoSet& set, const bool add) const
{
    // Note: wip. Select cells first
    boolList selectedCell(mesh_.nCells(), true);

    if (setName_.size() && setName_ != "none")
    {
        Info<< "    Loading subset " << setName_ << " to delimit search region."
            << endl;
        cellSet subSet(mesh_, setName_);

        selectedCell = false;
        forAllConstIter(labelHashSet, subSet, iter)
        {
            selectedCell[iter.key()] = true;
        }
    }


    unselectOutsideRegions(selectedCell);

    if (nErode_ > 0)
    {
        erode(selectedCell);
    }


    forAll(selectedCell, cellI)
    {
        if (selectedCell[cellI])
        {
            addOrDelete(set, cellI, add);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionsToCell::regionsToCell
(
    const polyMesh& mesh,
    const word& setName,
    const pointField& insidePoints,
    const label nErode
)
:
    topoSetSource(mesh),
    setName_(setName),
    insidePoints_(insidePoints),
    nErode_(nErode)
{}


Foam::regionsToCell::regionsToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    setName_(dict.lookupOrDefault<word>("set", "none")),
    insidePoints_
    (
        dict.found("insidePoints")
      ? dict.lookup("insidePoints")
      : dict.lookup("insidePoint")
    ),
    nErode_(dict.lookupOrDefault<label>("nErode", 0))
{}


Foam::regionsToCell::regionsToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    setName_(checkIs(is)),
    insidePoints_(checkIs(is)),
    nErode_(readLabel(checkIs(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionsToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (debug)
        {
            Info<< "    Adding all cells of connected region "
                << "containing points "
                << insidePoints_ << " ..." << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::REMOVE)
    {
        if (debug)
        {
            Info<< "    Removing all cells of connected region "
                << "containing points "
                << insidePoints_ << " ..." << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //
