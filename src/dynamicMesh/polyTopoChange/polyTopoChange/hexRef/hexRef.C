/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "hexRef.H"

#include "polyMesh.H"
#include "polyTopoChange.H"
#include "meshTools.H"
#include "polyAddFace.H"
#include "polyAddPoint.H"
#include "polyAddCell.H"
#include "polyModifyFace.H"
#include "syncTools.H"
#include "faceSet.H"
#include "cellSet.H"
#include "pointSet.H"
#include "OFstream.H"
#include "Time.H"
#include "FaceCellWave.H"
#include "mapDistributePolyMesh.H"
#include "refinementData.H"
#include "refinementDistanceData.H"
#include "degenerateMatcher.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hexRef, 0);
    defineRunTimeSelectionTable(hexRef, mesh);
    defineRunTimeSelectionTable(hexRef, levelsHist);
    defineRunTimeSelectionTable(hexRef, levels);

    //- Reduction class. If x and y are not equal assign value.
    template<label value>
    class ifEqEqOp
    {
        public:
        void operator()(label& x, const label y) const
        {
            x = (x == y) ? x : value;
        }
    };
}


// * * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * //

void Foam::hexRef::reorder
(
    const labelList& map,
    const label len,
    const label null,
    labelList& elems
)
{
    labelList newElems(len, null);

    forAll(elems, i)
    {
        label newI = map[i];

        if (newI >= len)
        {
            FatalErrorInFunction << abort(FatalError);
        }

        if (newI >= 0)
        {
            newElems[newI] = elems[i];
        }
    }

    elems.transfer(newElems);
}


void Foam::hexRef::getFaceInfo
(
    const label facei,
    label& patchID,
    label& zoneID,
    label& zoneFlip
) const
{
    patchID = -1;

    if (!mesh_.isInternalFace(facei))
    {
        patchID = mesh_.boundaryMesh().whichPatch(facei);
    }

    zoneID = mesh_.faceZones().whichZone(facei);

    zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh_.faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
    }
}


// Adds a face on top of existing facei.
Foam::label Foam::hexRef::addFace
(
    polyTopoChange& meshMod,
    const label facei,
    const face& newFace,
    const label own,
    const label nei
) const
{
    label patchID, zoneID, zoneFlip;

    getFaceInfo(facei, patchID, zoneID, zoneFlip);

    label newFacei = -1;

    if ((nei == -1) || (own < nei))
    {
        // Ordering ok.
        newFacei = meshMod.setAction
        (
            polyAddFace
            (
                newFace,                    // face
                own,                        // owner
                nei,                        // neighbour
                -1,                         // master point
                -1,                         // master edge
                facei,                      // master face for addition
                false,                      // flux flip
                patchID,                    // patch for face
                zoneID,                     // zone for face
                zoneFlip                    // face zone flip
            )
        );
    }
    else
    {
        // Reverse owner/neighbour
        newFacei = meshMod.setAction
        (
            polyAddFace
            (
                newFace.reverseFace(),      // face
                nei,                        // owner
                own,                        // neighbour
                -1,                         // master point
                -1,                         // master edge
                facei,                      // master face for addition
                false,                      // flux flip
                patchID,                    // patch for face
                zoneID,                     // zone for face
                zoneFlip                    // face zone flip
            )
        );
    }
    return newFacei;
}


// Adds an internal face from an edge. Assumes orientation correct.
// Problem is that the face is between four new vertices. So what do we provide
// as master? The only existing mesh item we have is the edge we have split.
// Have to be careful in only using it if it has internal faces since otherwise
// polyMeshMorph will complain (because it cannot generate a sensible mapping
// for the face)
Foam::label Foam::hexRef::addInternalFace
(
    polyTopoChange& meshMod,
    const label meshFacei,
    const label meshPointi,
    const face& newFace,
    const label own,
    const label nei
) const
{
/////////////////////////////////////////////////////////////////////////////////
// D. Deising & D. Rettenmaier
// Sets internalFaces falsly as masterFaces so internal face adressing gets mixed
// up. addInternalFace() is only called inside HexRef8, we expect no problems by
// commenting this case.

// TODO: Check if commenting this case breaks any other functionality

    // if (mesh_.isInternalFace(meshFacei))
    // {
    //    return meshMod.setAction
    //    (
    //        polyAddFace
    //        (
    //            newFace,                    // face
    //            own,                        // owner
    //            nei,                        // neighbour
    //            -1,                         // master point
    //            -1,                         // master edge
    //             meshFacei,                  // master face for addition
    //             false,                      // flux flip
    //             -1,                         // patch for face
    //             -1,                         // zone for face
    //             false                       // face zone flip
    //         )
    //     );
    // }
    // else
/////////////////////////////////////////////////////////////////////////////////
    {
        // Two choices:
        // - append (i.e. create out of nothing - will not be mapped)
        //   problem: field does not get mapped.
        // - inflate from point.
        //   problem: does interpolative mapping which constructs full
        //   volPointInterpolation!

        // For now create out of nothing

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


        ////- Inflate-from-point:
        //// Check if point has any internal faces we can use.
        //label masterPointi = -1;
        //
        //const labelList& pFaces = mesh_.pointFaces()[meshPointi];
        //
        //forAll(pFaces, i)
        //{
        //    if (mesh_.isInternalFace(pFaces[i]))
        //    {
        //        // meshPoint uses internal faces so ok to inflate from it
        //        masterPointi = meshPointi;
        //
        //        break;
        //    }
        //}
        //
        //return meshMod.setAction
        //(
        //    polyAddFace
        //    (
        //        newFace,                    // face
        //        own,                        // owner
        //        nei,                        // neighbour
        //        masterPointi,               // master point
        //        -1,                         // master edge
        //        -1,                         // master face for addition
        //        false,                      // flux flip
        //        -1,                         // patch for face
        //        -1,                         // zone for face
        //        false                       // face zone flip
        //    )
        //);
    }
}


// Modifies existing facei for either new owner/neighbour or new face points.
void Foam::hexRef::modFace
(
    polyTopoChange& meshMod,
    const label facei,
    const face& newFace,
    const label own,
    const label nei
) const
{
    label patchID, zoneID, zoneFlip;

    getFaceInfo(facei, patchID, zoneID, zoneFlip);

    if
    (
        (own != mesh_.faceOwner()[facei])
     || (
            mesh_.isInternalFace(facei)
         && (nei != mesh_.faceNeighbour()[facei])
        )
     || (newFace != mesh_.faces()[facei])
    )
    {
        if ((nei == -1) || (own < nei))
        {
            meshMod.setAction
            (
                polyModifyFace
                (
                    newFace,            // modified face
                    facei,              // label of face being modified
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
            meshMod.setAction
            (
                polyModifyFace
                (
                    newFace.reverseFace(),  // modified face
                    facei,                  // label of face being modified
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


// Bit complex way to determine the unrefined edge length.
Foam::scalar Foam::hexRef::getLevel0EdgeLength() const
{
    if (cellLevel_.size() != mesh_.nCells())
    {
        FatalErrorInFunction
            << "Number of cells in mesh:" << mesh_.nCells()
            << " does not equal size of cellLevel:" << cellLevel_.size()
            << endl
            << "This might be because of a restart with inconsistent cellLevel."
            << abort(FatalError);
    }

    // Determine minimum edge length per refinement level
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const scalar GREAT2 = sqr(GREAT);

    label nLevels = gMax(cellLevel_)+1;

    scalarField typEdgeLenSqr(nLevels, GREAT2);


    // 1. Look only at edges surrounded by cellLevel cells only.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        // Per edge the cellLevel of connected cells. -1 if not set,
        // labelMax if different levels, otherwise levels of connected cells.
        labelList edgeLevel(mesh_.nEdges(), -1);

        forAll(cellLevel_, celli)
        {
            const label cLevel = cellLevel_[celli];

            const labelList& cEdges = mesh_.cellEdges(celli);

            forAll(cEdges, i)
            {
                label edgeI = cEdges[i];

                if (edgeLevel[edgeI] == -1)
                {
                    edgeLevel[edgeI] = cLevel;
                }
                else if (edgeLevel[edgeI] == labelMax)
                {
                    // Already marked as on different cellLevels
                }
                else if (edgeLevel[edgeI] != cLevel)
                {
                    edgeLevel[edgeI] = labelMax;
                }
            }
        }

        // Make sure that edges with different levels on different processors
        // are also marked. Do the same test (edgeLevel != cLevel) on coupled
        // edges.
        syncTools::syncEdgeList
        (
            mesh_,
            edgeLevel,
            ifEqEqOp<labelMax>(),
            labelMin
        );

        // Now use the edgeLevel with a valid value to determine the
        // length per level.
        forAll(edgeLevel, edgeI)
        {
            const label eLevel = edgeLevel[edgeI];

            if (eLevel >= 0 && eLevel < labelMax)
            {
                const edge& e = mesh_.edges()[edgeI];

                scalar edgeLenSqr = magSqr(e.vec(mesh_.points()));

                typEdgeLenSqr[eLevel] = min(typEdgeLenSqr[eLevel], edgeLenSqr);
            }
        }
    }

    // Get the minimum per level over all processors. Note minimum so if
    // cells are not cubic we use the smallest edge side.
    Pstream::listCombineGather(typEdgeLenSqr, minEqOp<scalar>());
    Pstream::listCombineScatter(typEdgeLenSqr);

    if (debug)
    {
        Pout<< "hexRef::getLevel0EdgeLength() :"
            << " After phase1: Edgelengths (squared) per refinementlevel:"
            << typEdgeLenSqr << endl;
    }


    // 2. For any levels where we haven't determined a valid length yet
    //    use any surrounding cell level. Here we use the max so we don't
    //    pick up levels between celllevel and higher celllevel (will have
    //    edges sized according to highest celllevel)
    //    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    scalarField maxEdgeLenSqr(nLevels, -GREAT2);

    forAll(cellLevel_, celli)
    {
        const label cLevel = cellLevel_[celli];

        const labelList& cEdges = mesh_.cellEdges(celli);

        forAll(cEdges, i)
        {
            const edge& e = mesh_.edges()[cEdges[i]];

            scalar edgeLenSqr = magSqr(e.vec(mesh_.points()));

            maxEdgeLenSqr[cLevel] = max(maxEdgeLenSqr[cLevel], edgeLenSqr);
        }
    }

    Pstream::listCombineGather(maxEdgeLenSqr, maxEqOp<scalar>());
    Pstream::listCombineScatter(maxEdgeLenSqr);

    if (debug)
    {
        Pout<< "hexRef::getLevel0EdgeLength() :"
            << " Crappy Edgelengths (squared) per refinementlevel:"
            << maxEdgeLenSqr << endl;
    }


    // 3. Combine the two sets of lengths
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(typEdgeLenSqr, levelI)
    {
        if (typEdgeLenSqr[levelI] == GREAT2 && maxEdgeLenSqr[levelI] >= 0)
        {
            typEdgeLenSqr[levelI] = maxEdgeLenSqr[levelI];
        }
    }

    if (debug)
    {
        Pout<< "hexRef::getLevel0EdgeLength() :"
            << " Final Edgelengths (squared) per refinementlevel:"
            << typEdgeLenSqr << endl;
    }

    // Find lowest level present
    scalar level0Size = -1;

    forAll(typEdgeLenSqr, levelI)
    {
        scalar lenSqr = typEdgeLenSqr[levelI];

        if (lenSqr < GREAT2)
        {
            level0Size = Foam::sqrt(lenSqr)*(1<<levelI);

            if (debug)
            {
                Pout<< "hexRef::getLevel0EdgeLength() :"
                    << " For level:" << levelI
                    << " have edgeLen:" << Foam::sqrt(lenSqr)
                    << " with equivalent level0 len:" << level0Size
                    << endl;
            }
            break;
        }
    }

    if (level0Size == -1)
    {
        FatalErrorInFunction
            << "Problem : typEdgeLenSqr:" << typEdgeLenSqr << abort(FatalError);
    }

    return level0Size;
}


// Check whether pointi is an anchor on celli.
// If it is not check whether any other point on the face is an anchor cell.
Foam::label Foam::hexRef::getAnchorCell
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const label celli,
    const label facei,
    const label pointi
) const
{
    if (cellAnchorPoints[celli].size())
    {
        label index = findIndex(cellAnchorPoints[celli], pointi);

        if (index != -1)
        {
            return cellAddedCells[celli][index];
        }


        // pointi is not an anchor cell.
        // Maybe we are already a refined face so check all the face
        // vertices.
        const face& f = mesh_.faces()[facei];

        forAll(f, fp)
        {
            label index = findIndex(cellAnchorPoints[celli], f[fp]);

            if (index != -1)
            {
                return cellAddedCells[celli][index];
            }
        }

        // Problem.
        dumpCell(celli);
        Perr<< "cell:" << celli << " anchorPoints:" << cellAnchorPoints[celli]
            << endl;

        FatalErrorInFunction
            << "Could not find point " << pointi
            << " in the anchorPoints for cell " << celli << endl
            << "Does your original mesh obey the 2:1 constraint and"
            << " did you use consistentRefinement to make your cells to refine"
            << " obey this constraint as well?"
            << abort(FatalError);

        return -1;
    }
    else
    {
        return celli;
    }
}


// Get new owner and neighbour
void Foam::hexRef::getFaceNeighbours
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const label facei,
    const label pointi,

    label& own,
    label& nei
) const
{
    // Is owner split?
    own = getAnchorCell
    (
        cellAnchorPoints,
        cellAddedCells,
        mesh_.faceOwner()[facei],
        facei,
        pointi
    );

    if (mesh_.isInternalFace(facei))
    {
        nei = getAnchorCell
        (
            cellAnchorPoints,
            cellAddedCells,
            mesh_.faceNeighbour()[facei],
            facei,
            pointi
        );
    }
    else
    {
        nei = -1;
    }
}


// Get point with the lowest pointLevel
Foam::label Foam::hexRef::findMinLevel(const labelList& f) const
{
    label minLevel = labelMax;
    label minFp = -1;

    forAll(f, fp)
    {
        label level = pointLevel_[f[fp]];

        if (level < minLevel)
        {
            minLevel = level;
            minFp = fp;
        }
    }

    return minFp;
}


// Get point with the highest pointLevel
Foam::label Foam::hexRef::findMaxLevel(const labelList& f) const
{
    label maxLevel = labelMin;
    label maxFp = -1;

    forAll(f, fp)
    {
        label level = pointLevel_[f[fp]];

        if (level > maxLevel)
        {
            maxLevel = level;
            maxFp = fp;
        }
    }

    return maxFp;
}


Foam::label Foam::hexRef::countAnchors
(
    const labelList& f,
    const label anchorLevel
) const
{
    label nAnchors = 0;

    forAll(f, fp)
    {
        if (pointLevel_[f[fp]] <= anchorLevel)
        {
            nAnchors++;
        }
    }
    return nAnchors;
}


void Foam::hexRef::dumpCell(const label celli) const
{
    OFstream str(mesh_.time().path()/"cell_" + Foam::name(celli) + ".obj");
    Pout<< "hexRef : Dumping cell as obj to " << str.name() << endl;

    const cell& cFaces = mesh_.cells()[celli];

    Map<label> pointToObjVert;
    label objVertI = 0;

    forAll(cFaces, i)
    {
        const face& f = mesh_.faces()[cFaces[i]];

        forAll(f, fp)
        {
            if (pointToObjVert.insert(f[fp], objVertI))
            {
                meshTools::writeOBJ(str, mesh_.points()[f[fp]]);
                objVertI++;
            }
        }
    }

    forAll(cFaces, i)
    {
        const face& f = mesh_.faces()[cFaces[i]];

        forAll(f, fp)
        {
            label pointi = f[fp];
            label nexPointi = f[f.fcIndex(fp)];

            str << "l " << pointToObjVert[pointi]+1
                << ' ' << pointToObjVert[nexPointi]+1 << nl;
        }
    }
}


// Find point with certain pointLevel. Skip any higher levels.
Foam::label Foam::hexRef::findLevel
(
    const label facei,
    const face& f,
    const label startFp,
    const bool searchForward,
    const label wantedLevel
) const
{
    label fp = startFp;

    forAll(f, i)
    {
        label pointi = f[fp];

        if (pointLevel_[pointi] < wantedLevel)
        {
            dumpCell(mesh_.faceOwner()[facei]);
            if (mesh_.isInternalFace(facei))
            {
                dumpCell(mesh_.faceNeighbour()[facei]);
            }

            FatalErrorInFunction
                << "face:" << f
                << " level:" << UIndirectList<label>(pointLevel_, f)()
                << " startFp:" << startFp
                << " wantedLevel:" << wantedLevel
                << abort(FatalError);
        }
        else if (pointLevel_[pointi] == wantedLevel)
        {
            return fp;
        }

        if (searchForward)
        {
            fp = f.fcIndex(fp);
        }
        else
        {
            fp = f.rcIndex(fp);
        }
    }

    dumpCell(mesh_.faceOwner()[facei]);
    if (mesh_.isInternalFace(facei))
    {
        dumpCell(mesh_.faceNeighbour()[facei]);
    }

    FatalErrorInFunction
        << "face:" << f
        << " level:" << UIndirectList<label>(pointLevel_, f)()
        << " startFp:" << startFp
        << " wantedLevel:" << wantedLevel
        << abort(FatalError);

    return -1;
}


// Gets cell level such that the face has four points <= level.
Foam::label Foam::hexRef::faceLevel(const label facei) const
{
    const face& f = mesh_.faces()[facei];

    if (f.size() <= 4)
    {
        return pointLevel_[f[findMaxLevel(f)]];
    }
    else
    {
        label ownLevel = cellLevel_[mesh_.faceOwner()[facei]];

        if (countAnchors(f, ownLevel) == 4)
        {
            return ownLevel;
        }
        else if (countAnchors(f, ownLevel+1) == 4)
        {
            return ownLevel+1;
        }
        else
        {
            return -1;
        }
    }
}


void Foam::hexRef::checkInternalOrientation
(
    polyTopoChange& meshMod,
    const label celli,
    const label facei,
    const point& ownPt,
    const point& neiPt,
    const face& newFace
)
{
    face compactFace(identity(newFace.size()));
    pointField compactPoints(meshMod.points(), newFace);

    vector n(compactFace.normal(compactPoints));

    vector dir(neiPt - ownPt);

    if ((dir & n) < 0)
    {
        FatalErrorInFunction
            << "cell:" << celli << " old face:" << facei
            << " newFace:" << newFace << endl
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " neiPt:" << neiPt
            << abort(FatalError);
    }

    vector fcToOwn(compactFace.centre(compactPoints) - ownPt);

    scalar s = (fcToOwn&n) / (dir&n);

    if (s < 0.1 || s > 0.9)
    {
        FatalErrorInFunction
            << "cell:" << celli << " old face:" << facei
            << " newFace:" << newFace << endl
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " neiPt:" << neiPt
            << " s:" << s
            << abort(FatalError);
    }
}


void Foam::hexRef::checkBoundaryOrientation
(
    polyTopoChange& meshMod,
    const label celli,
    const label facei,
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
            << "cell:" << celli << " old face:" << facei
            << " newFace:" << newFace
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " boundaryPt:" << boundaryPt
            << abort(FatalError);
    }

    vector fcToOwn(compactFace.centre(compactPoints) - ownPt);

    scalar s = (fcToOwn&dir) / magSqr(dir);

    if (s < 0.7 || s > 1.3)
    {
        WarningInFunction
            << "cell:" << celli << " old face:" << facei
            << " newFace:" << newFace
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " boundaryPt:" << boundaryPt
            << " s:" << s
            << endl;
            //<< abort(FatalError);
    }
}


// If p0 and p1 are existing vertices check if edge is split and insert
// splitPoint.
void Foam::hexRef::insertEdgeSplit
(
    const labelList& edgeMidPoint,
    const label p0,
    const label p1,
    DynamicList<label>& verts
) const
{
    if (p0 < mesh_.nPoints() && p1 < mesh_.nPoints())
    {
        label edgeI = meshTools::findEdge(mesh_, p0, p1);

        if (edgeI != -1 && edgeMidPoint[edgeI] != -1)
        {
            verts.append(edgeMidPoint[edgeI]);
        }
    }
}


void Foam::hexRef::walkFaceToMid
(
    const labelList& edgeMidPoint,
    const label cLevel,
    const label facei,
    const label startFp,
    DynamicList<label>& faceVerts
) const
{
    const face& f = mesh_.faces()[facei];
    const labelList& fEdges = mesh_.faceEdges(facei);

    label fp = startFp;

    // Starting from fp store all (1 or 2) vertices until where the face
    // gets split

    while (true)
    {
        if (edgeMidPoint[fEdges[fp]] >= 0)
        {
            faceVerts.append(edgeMidPoint[fEdges[fp]]);
        }

        fp = f.fcIndex(fp);

        if (pointLevel_[f[fp]] <= cLevel)
        {
            // Next anchor. Have already append split point on edge in code
            // above.
            return;
        }
        else if (pointLevel_[f[fp]] == cLevel+1)
        {
            // Mid level
            faceVerts.append(f[fp]);

            return;
        }
        else if (pointLevel_[f[fp]] == cLevel+2)
        {
            // Store and continue to cLevel+1.
            faceVerts.append(f[fp]);
        }
    }
}


// Same as walkFaceToMid but now walk back.
void Foam::hexRef::walkFaceFromMid
(
    const labelList& edgeMidPoint,
    const label cLevel,
    const label facei,
    const label startFp,
    DynamicList<label>& faceVerts
) const
{
    const face& f = mesh_.faces()[facei];
    const labelList& fEdges = mesh_.faceEdges(facei);

    label fp = f.rcIndex(startFp);

    while (true)
    {
        if (pointLevel_[f[fp]] <= cLevel)
        {
            // anchor.
            break;
        }
        else if (pointLevel_[f[fp]] == cLevel+1)
        {
            // Mid level
            faceVerts.append(f[fp]);
            break;
        }
        else if (pointLevel_[f[fp]] == cLevel+2)
        {
            // Continue to cLevel+1.
        }
        fp = f.rcIndex(fp);
    }

    // Store
    while (true)
    {
        if (edgeMidPoint[fEdges[fp]] >= 0)
        {
            faceVerts.append(edgeMidPoint[fEdges[fp]]);
        }

        fp = f.fcIndex(fp);

        if (fp == startFp)
        {
            break;
        }
        faceVerts.append(f[fp]);
    }
}


// Updates refineCell (cells marked for refinement) so across all faces
// there will be 2:1 consistency after refinement.
Foam::label Foam::hexRef::faceConsistentRefinement
(
    const bool maxSet,
    PackedBoolList& refineCell
) const
{
    label nChanged = 0;

    // Internal faces.
    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        label own = mesh_.faceOwner()[facei];
        label ownLevel = cellLevel_[own] + refineCell.get(own);

        label nei = mesh_.faceNeighbour()[facei];
        label neiLevel = cellLevel_[nei] + refineCell.get(nei);

        if (ownLevel > (neiLevel+1))
        {
            if (maxSet)
            {
                refineCell.set(nei);
            }
            else
            {
                refineCell.unset(own);
            }
            nChanged++;
        }
        else if (neiLevel > (ownLevel+1))
        {
            if (maxSet)
            {
                refineCell.set(own);
            }
            else
            {
                refineCell.unset(nei);
            }
            nChanged++;
        }
    }


    // Coupled faces. Swap owner level to get neighbouring cell level.
    // (only boundary faces of neiLevel used)
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

    forAll(neiLevel, i)
    {
        label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];

        neiLevel[i] = cellLevel_[own] + refineCell.get(own);
    }

    // Swap to neighbour
    syncTools::swapBoundaryFaceList(mesh_, neiLevel);

    // Now we have neighbour value see which cells need refinement
    forAll(neiLevel, i)
    {
        label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];
        label ownLevel = cellLevel_[own] + refineCell.get(own);

        if (ownLevel > (neiLevel[i]+1))
        {
            if (!maxSet)
            {
                refineCell.unset(own);
                nChanged++;
            }
        }
        else if (neiLevel[i] > (ownLevel+1))
        {
            if (maxSet)
            {
                refineCell.set(own);
                nChanged++;
            }
        }
    }

    return nChanged;
}


// Debug: check if wanted refinement is compatible with 2:1
void Foam::hexRef::checkWantedRefinementLevels
(
    const labelList& cellsToRefine
) const
{
    PackedBoolList refineCell(mesh_.nCells());
    forAll(cellsToRefine, i)
    {
        refineCell.set(cellsToRefine[i]);
    }

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        label own = mesh_.faceOwner()[facei];
        label ownLevel = cellLevel_[own] + refineCell.get(own);

        label nei = mesh_.faceNeighbour()[facei];
        label neiLevel = cellLevel_[nei] + refineCell.get(nei);

        if (mag(ownLevel-neiLevel) > 1)
        {
            dumpCell(own);
            dumpCell(nei);
            FatalErrorInFunction
                << "cell:" << own
                << " current level:" << cellLevel_[own]
                << " level after refinement:" << ownLevel
                << nl
                << "neighbour cell:" << nei
                << " current level:" << cellLevel_[nei]
                << " level after refinement:" << neiLevel
                << nl
                << "which does not satisfy 2:1 constraints anymore."
                << abort(FatalError);
        }
    }

    // Coupled faces. Swap owner level to get neighbouring cell level.
    // (only boundary faces of neiLevel used)
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

    forAll(neiLevel, i)
    {
        label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];

        neiLevel[i] = cellLevel_[own] + refineCell.get(own);
    }

    // Swap to neighbour
    syncTools::swapBoundaryFaceList(mesh_, neiLevel);

    // Now we have neighbour value see which cells need refinement
    forAll(neiLevel, i)
    {
        label facei = i + mesh_.nInternalFaces();

        label own = mesh_.faceOwner()[facei];
        label ownLevel = cellLevel_[own] + refineCell.get(own);

        if (mag(ownLevel - neiLevel[i]) > 1)
        {
            label patchi = mesh_.boundaryMesh().whichPatch(facei);

            dumpCell(own);
            FatalErrorInFunction
                << "Celllevel does not satisfy 2:1 constraint."
                << " On coupled face "
                << facei
                << " on patch " << patchi << " "
                << mesh_.boundaryMesh()[patchi].name()
                << " owner cell " << own
                << " current level:" << cellLevel_[own]
                << " level after refinement:" << ownLevel
                << nl
                << " (coupled) neighbour cell will get refinement "
                << neiLevel[i]
                << abort(FatalError);
        }
    }
}


// Set instance for mesh files
void Foam::hexRef::setInstance(const fileName& inst)
{
    if (debug)
    {
        Pout<< "hexRef::setInstance(const fileName& inst) : "
            << "Resetting file instance to " << inst << endl;
    }

    cellLevel_.instance() = inst;
    pointLevel_.instance() = inst;
    level0Edge_.instance() = inst;
    history_.instance() = inst;
}


void Foam::hexRef::collectLevelPoints
(
    const labelList& f,
    const label level,
    DynamicList<label>& points
) const
{
    forAll(f, fp)
    {
        if (pointLevel_[f[fp]] <= level)
        {
            points.append(f[fp]);
        }
    }
}


void Foam::hexRef::collectLevelPoints
(
    const labelList& meshPoints,
    const labelList& f,
    const label level,
    DynamicList<label>& points
) const
{
    forAll(f, fp)
    {
        label pointi = meshPoints[f[fp]];
        if (pointLevel_[pointi] <= level)
        {
            points.append(pointi);
        }
    }
}


// Return true if we've found 6 quads. faces guaranteed to be outwards pointing.
bool Foam::hexRef::matchHexShape
(
    const label celli,
    const label cellLevel,
    DynamicList<face>& quads
) const
{
    const cell& cFaces = mesh_.cells()[celli];

    // Work arrays
    DynamicList<label> verts(4);
    quads.clear();


    // 1. pick up any faces with four cellLevel points

    forAll(cFaces, i)
    {
        label facei = cFaces[i];
        const face& f = mesh_.faces()[facei];

        verts.clear();
        collectLevelPoints(f, cellLevel, verts);
        if (verts.size() == 4)
        {
            if (mesh_.faceOwner()[facei] != celli)
            {
                reverse(verts);
            }
            quads.append(face(0));
            labelList& quadVerts = quads.last();
            quadVerts.transfer(verts);
        }
    }


    if (quads.size() < 6)
    {
        Map<labelList> pointFaces(2*cFaces.size());

        forAll(cFaces, i)
        {
            label facei = cFaces[i];
            const face& f = mesh_.faces()[facei];

            // Pick up any faces with only one level point.
            // See if there are four of these where the commont point
            // is a level+1 point. This common point is then the mid of
            // a split face.

            verts.clear();
            collectLevelPoints(f, cellLevel, verts);
            if (verts.size() == 1)
            {
                // Add to pointFaces for any level+1 point (this might be
                // a midpoint of a split face)
                forAll(f, fp)
                {
                    label pointi = f[fp];
                    if (pointLevel_[pointi] == cellLevel+1)
                    {
                        Map<labelList>::iterator iter =
                            pointFaces.find(pointi);
                        if (iter != pointFaces.end())
                        {
                            labelList& pFaces = iter();
                            if (findIndex(pFaces, facei) == -1)
                            {
                                pFaces.append(facei);
                            }
                        }
                        else
                        {
                            pointFaces.insert
                            (
                                pointi,
                                labelList(1, facei)
                            );
                        }
                    }
                }
            }
        }

        // 2. Check if we've collected any midPoints.
        forAllConstIter(Map<labelList>, pointFaces, iter)
        {
            const labelList& pFaces = iter();

            if (pFaces.size() == 4)
            {
                // Collect and orient.
                faceList fourFaces(pFaces.size());
                forAll(pFaces, pFacei)
                {
                    label facei = pFaces[pFacei];
                    const face& f = mesh_.faces()[facei];
                    if (mesh_.faceOwner()[facei] == celli)
                    {
                        fourFaces[pFacei] = f;
                    }
                    else
                    {
                        fourFaces[pFacei] = f.reverseFace();
                    }
                }

                primitivePatch bigFace
                (
                    SubList<face>(fourFaces, fourFaces.size()),
                    mesh_.points()
                );
                const labelListList& edgeLoops = bigFace.edgeLoops();

                if (edgeLoops.size() == 1)
                {
                    // Collect the 4 cellLevel points
                    verts.clear();
                    collectLevelPoints
                    (
                        bigFace.meshPoints(),
                        bigFace.edgeLoops()[0],
                        cellLevel,
                        verts
                    );

                    if (verts.size() == 4)
                    {
                        quads.append(face(0));
                        labelList& quadVerts = quads.last();
                        quadVerts.transfer(verts);
                    }
                }
            }
        }
    }

    return (quads.size() == 6);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh, read refinement data
Foam::hexRef::hexRef(const polyMesh& mesh, const bool readHistory)
:
    mesh_(mesh),
    cellLevel_
    (
        IOobject
        (
            "cellLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        labelList(mesh_.nCells(), 0)
    ),
    pointLevel_
    (
        IOobject
        (
            "pointLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        labelList(mesh_.nPoints(), 0)
    ),
    level0Edge_
    (
        IOobject
        (
            "level0Edge",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar("level0Edge", dimLength, getLevel0EdgeLength())
    ),
    history_
    (
        IOobject
        (
            "refinementHistory",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        // All cells visible if not read or readHistory = false
        (readHistory ? mesh_.nCells() : 0)
    ),
    faceRemover_(mesh_, GREAT),     // merge boundary faces wherever possible
    savedPointLevel_(0),
    savedCellLevel_(0)
{
    if (readHistory)
    {
        // Make sure we don't use the master-only reading. Bit of a hack for
        // now.
        regIOobject::fileCheckTypes oldType =
            regIOobject::fileModificationChecking;
        regIOobject::fileModificationChecking = regIOobject::timeStamp;
        history_.readOpt() = IOobject::READ_IF_PRESENT;
        if (history_.headerOk())
        {
            history_.read();
        }
        regIOobject::fileModificationChecking = oldType;
    }

    if (history_.active() && history_.visibleCells().size() != mesh_.nCells())
    {
        FatalErrorInFunction
            << "History enabled but number of visible cells "
            << history_.visibleCells().size() << " in "
            << history_.objectPath()
            << " is not equal to the number of cells in the mesh "
            << mesh_.nCells()
            << abort(FatalError);
    }

    if
    (
        cellLevel_.size() != mesh_.nCells()
     || pointLevel_.size() != mesh_.nPoints()
    )
    {
        FatalErrorInFunction
            << "Restarted from inconsistent cellLevel or pointLevel files."
            << endl
            << "cellLevel file " << cellLevel_.objectPath() << endl
            << "pointLevel file " << pointLevel_.objectPath() << endl
            << "Number of cells in mesh:" << mesh_.nCells()
            << " does not equal size of cellLevel:" << cellLevel_.size() << endl
            << "Number of points in mesh:" << mesh_.nPoints()
            << " does not equal size of pointLevel:" << pointLevel_.size()
            << abort(FatalError);
    }


    // Check refinement levels for consistency
    checkRefinementLevels(-1, labelList(0));


    // Check initial mesh for consistency

    //if (debug)
    {
        checkMesh();
    }
}


// Construct from components
Foam::hexRef::hexRef
(
    const polyMesh& mesh,
    const labelList& cellLevel,
    const labelList& pointLevel,
    const refinementHistory& history,
    const scalar level0Edge
)
:
    mesh_(mesh),
    cellLevel_
    (
        IOobject
        (
            "cellLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cellLevel
    ),
    pointLevel_
    (
        IOobject
        (
            "pointLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pointLevel
    ),
    level0Edge_
    (
        IOobject
        (
            "level0Edge",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dimensionedScalar
        (
            "level0Edge",
            dimLength,
            (level0Edge >= 0 ? level0Edge : getLevel0EdgeLength())
        )
    ),
    history_
    (
        IOobject
        (
            "refinementHistory",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        history
    ),
    faceRemover_(mesh_, GREAT),     // merge boundary faces wherever possible
    savedPointLevel_(0),
    savedCellLevel_(0)
{
    if (history_.active() && history_.visibleCells().size() != mesh_.nCells())
    {
        FatalErrorInFunction
            << "History enabled but number of visible cells in it "
            << history_.visibleCells().size()
            << " is not equal to the number of cells in the mesh "
            << mesh_.nCells() << abort(FatalError);
    }

    if
    (
        cellLevel_.size() != mesh_.nCells()
     || pointLevel_.size() != mesh_.nPoints()
    )
    {
        FatalErrorInFunction
            << "Incorrect cellLevel or pointLevel size." << endl
            << "Number of cells in mesh:" << mesh_.nCells()
            << " does not equal size of cellLevel:" << cellLevel_.size() << endl
            << "Number of points in mesh:" << mesh_.nPoints()
            << " does not equal size of pointLevel:" << pointLevel_.size()
            << abort(FatalError);
    }

    // Check refinement levels for consistency
    checkRefinementLevels(-1, labelList(0));


    // Check initial mesh for consistency

    //if (debug)
    {
        checkMesh();
    }
}


// Construct from components
Foam::hexRef::hexRef
(
    const polyMesh& mesh,
    const labelList& cellLevel,
    const labelList& pointLevel,
    const scalar level0Edge
)
:
    mesh_(mesh),
    cellLevel_
    (
        IOobject
        (
            "cellLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cellLevel
    ),
    pointLevel_
    (
        IOobject
        (
            "pointLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pointLevel
    ),
    level0Edge_
    (
        IOobject
        (
            "level0Edge",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dimensionedScalar
        (
            "level0Edge",
            dimLength,
            (level0Edge >= 0 ? level0Edge : getLevel0EdgeLength())
        )
    ),
    history_
    (
        IOobject
        (
            "refinementHistory",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        List<refinementHistory::splitCell8>(0),
        labelList(0),
        false
    ),
    faceRemover_(mesh_, GREAT),     // merge boundary faces wherever possible
    savedPointLevel_(0),
    savedCellLevel_(0)
{
    if
    (
        cellLevel_.size() != mesh_.nCells()
     || pointLevel_.size() != mesh_.nPoints()
    )
    {
        FatalErrorInFunction
            << "Incorrect cellLevel or pointLevel size." << endl
            << "Number of cells in mesh:" << mesh_.nCells()
            << " does not equal size of cellLevel:" << cellLevel_.size() << endl
            << "Number of points in mesh:" << mesh_.nPoints()
            << " does not equal size of pointLevel:" << pointLevel_.size()
            << abort(FatalError);
    }

    // Check refinement levels for consistency
    checkRefinementLevels(-1, labelList(0));

    // Check initial mesh for consistency

    //if (debug)
    {
        checkMesh();
    }
}

// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::hexRef::~hexRef()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::hexRef::consistentRefinement
(
    const labelList& cellsToRefine,
    const bool maxSet
) const
{
    // Loop, modifying cellsToRefine, until no more changes to due to 2:1
    // conflicts.
    // maxSet = false : unselect cells to refine
    // maxSet = true  : select cells to refine

    // Go to straight boolList.
    PackedBoolList refineCell(mesh_.nCells());
    forAll(cellsToRefine, i)
    {
        refineCell.set(cellsToRefine[i]);
    }

    while (true)
    {
        label nChanged = faceConsistentRefinement(maxSet, refineCell);

        reduce(nChanged, sumOp<label>());

        if (debug)
        {
            Pout<< "hexRef::consistentRefinement : Changed " << nChanged
                << " refinement levels due to 2:1 conflicts."
                << endl;
        }

        if (nChanged == 0)
        {
            break;
        }
    }


    // Convert back to labelList.
    label nRefined = 0;

    forAll(refineCell, celli)
    {
        if (refineCell.get(celli))
        {
            nRefined++;
        }
    }

    labelList newCellsToRefine(nRefined);
    nRefined = 0;

    forAll(refineCell, celli)
    {
        if (refineCell.get(celli))
        {
            newCellsToRefine[nRefined++] = celli;
        }
    }

    if (debug)
    {
        checkWantedRefinementLevels(newCellsToRefine);
    }

    return newCellsToRefine;
}


// Given a list of cells to refine determine additional cells to refine
// such that the overall refinement:
// - satisfies maxFaceDiff (e.g. 2:1) across neighbouring faces
// - satisfies maxPointDiff (e.g. 4:1) across selected point connected
//   cells. This is used to ensure that e.g. cells on the surface are not
//   point connected to cells which are 8 times smaller.
Foam::labelList Foam::hexRef::consistentSlowRefinement
(
    const label maxFaceDiff,
    const labelList& cellsToRefine,
    const labelList& facesToCheck,
    const label maxPointDiff,
    const labelList& pointsToCheck
) const
{
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();


    if (maxFaceDiff <= 0)
    {
        FatalErrorInFunction
            << "Illegal maxFaceDiff " << maxFaceDiff << nl
            << "Value should be >= 1" << exit(FatalError);
    }


    // Bit tricky. Say we want a distance of three cells between two
    // consecutive refinement levels. This is done by using FaceCellWave to
    // transport out the new refinement level. It gets decremented by one
    // every cell it crosses so if we initialize it to maxFaceDiff
    // we will get a field everywhere that tells us whether an unselected cell
    // needs refining as well.


    // Initial information about (distance to) cellLevel on all cells
    List<refinementData> allCellInfo(mesh_.nCells());

    // Initial information about (distance to) cellLevel on all faces
    List<refinementData> allFaceInfo(mesh_.nFaces());

    forAll(allCellInfo, celli)
    {
        // maxFaceDiff since refinementData counts both
        // faces and cells.
        allCellInfo[celli] = refinementData
        (
            maxFaceDiff*(cellLevel_[celli]+1),// when cell is to be refined
            maxFaceDiff*cellLevel_[celli]     // current level
        );
    }

    // Cells to be refined will have cellLevel+1
    forAll(cellsToRefine, i)
    {
        label celli = cellsToRefine[i];

        allCellInfo[celli].count() = allCellInfo[celli].refinementCount();
    }


    // Labels of seed faces
    DynamicList<label> seedFaces(mesh_.nFaces()/100);
    // refinementLevel data on seed faces
    DynamicList<refinementData> seedFacesInfo(mesh_.nFaces()/100);

    // Dummy additional info for FaceCellWave
    int dummyTrackData = 0;


    // Additional buffer layer thickness by changing initial count. Usually
    // this happens on boundary faces. Bit tricky. Use allFaceInfo to mark
    // off thus marked faces so they're skipped in the next loop.
    forAll(facesToCheck, i)
    {
        label facei = facesToCheck[i];

        if (allFaceInfo[facei].valid(dummyTrackData))
        {
            // Can only occur if face has already gone through loop below.
            FatalErrorInFunction
                << "Argument facesToCheck seems to have duplicate entries!"
                << endl
                << "face:" << facei << " occurs at positions "
                << findIndices(facesToCheck, facei)
                << abort(FatalError);
        }


        const refinementData& ownData = allCellInfo[faceOwner[facei]];

        if (mesh_.isInternalFace(facei))
        {
            // Seed face if neighbouring cell (after possible refinement)
            // will be refined one more than the current owner or neighbour.

            const refinementData& neiData = allCellInfo[faceNeighbour[facei]];

            label faceCount;
            label faceRefineCount;
            if (neiData.count() > ownData.count())
            {
                faceCount = neiData.count() + maxFaceDiff;
                faceRefineCount = faceCount + maxFaceDiff;
            }
            else
            {
                faceCount = ownData.count() + maxFaceDiff;
                faceRefineCount = faceCount + maxFaceDiff;
            }

            seedFaces.append(facei);
            seedFacesInfo.append
            (
                refinementData
                (
                    faceRefineCount,
                    faceCount
                )
            );
            allFaceInfo[facei] = seedFacesInfo.last();
        }
        else
        {
            label faceCount = ownData.count() + maxFaceDiff;
            label faceRefineCount = faceCount + maxFaceDiff;

            seedFaces.append(facei);
            seedFacesInfo.append
            (
                refinementData
                (
                    faceRefineCount,
                    faceCount
                )
            );
            allFaceInfo[facei] = seedFacesInfo.last();
        }
    }


    // Just seed with all faces inbetween different refinement levels for now
    // (alternatively only seed faces on cellsToRefine but that gives problems
    //  if no cells to refine)
    forAll(faceNeighbour, facei)
    {
        // Check if face already handled in loop above
        if (!allFaceInfo[facei].valid(dummyTrackData))
        {
            label own = faceOwner[facei];
            label nei = faceNeighbour[facei];

            // Seed face with transported data from highest cell.

            if (allCellInfo[own].count() > allCellInfo[nei].count())
            {
                allFaceInfo[facei].updateFace
                (
                    mesh_,
                    facei,
                    own,
                    allCellInfo[own],
                    FaceCellWave<refinementData, int>::propagationTol(),
                    dummyTrackData
                );
                seedFaces.append(facei);
                seedFacesInfo.append(allFaceInfo[facei]);
            }
            else if (allCellInfo[own].count() < allCellInfo[nei].count())
            {
                allFaceInfo[facei].updateFace
                (
                    mesh_,
                    facei,
                    nei,
                    allCellInfo[nei],
                    FaceCellWave<refinementData, int>::propagationTol(),
                    dummyTrackData
                );
                seedFaces.append(facei);
                seedFacesInfo.append(allFaceInfo[facei]);
            }
        }
    }

    // Seed all boundary faces with owner value. This is to make sure that
    // they are visited (probably only important for coupled faces since
    // these need to be visited from both sides)
    for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); facei++)
    {
        // Check if face already handled in loop above
        if (!allFaceInfo[facei].valid(dummyTrackData))
        {
            label own = faceOwner[facei];

            // Seed face with transported data from owner.
            refinementData faceData;
            faceData.updateFace
            (
                mesh_,
                facei,
                own,
                allCellInfo[own],
                FaceCellWave<refinementData, int>::propagationTol(),
                dummyTrackData
            );
            seedFaces.append(facei);
            seedFacesInfo.append(faceData);
        }
    }


    // face-cell-face transport engine
    FaceCellWave<refinementData, int> levelCalc
    (
        mesh_,
        allFaceInfo,
        allCellInfo,
        dummyTrackData
    );

    while (true)
    {
        if (debug)
        {
            Pout<< "hexRef::consistentSlowRefinement : Seeded "
                << seedFaces.size() << " faces between cells with different"
                << " refinement level." << endl;
        }

        // Set seed faces
        levelCalc.setFaceInfo(seedFaces.shrink(), seedFacesInfo.shrink());
        seedFaces.clear();
        seedFacesInfo.clear();

        // Iterate until no change. Now 2:1 face difference should be satisfied
        levelCalc.iterate(mesh_.globalData().nTotalFaces()+1);


        // Now check point-connected cells (face-connected cells already ok):
        // - get per point max of connected cells
        // - sync across coupled points
        // - check cells against above point max

        if (maxPointDiff == -1)
        {
            // No need to do any point checking.
            break;
        }

        // Determine per point the max cell level. (done as count, not
        // as cell level purely for ease)
        labelList maxPointCount(mesh_.nPoints(), 0);

        forAll(maxPointCount, pointi)
        {
            label& pLevel = maxPointCount[pointi];

            const labelList& pCells = mesh_.pointCells(pointi);

            forAll(pCells, i)
            {
                pLevel = max(pLevel, allCellInfo[pCells[i]].count());
            }
        }

        // Sync maxPointCount to neighbour
        syncTools::syncPointList
        (
            mesh_,
            maxPointCount,
            maxEqOp<label>(),
            labelMin            // null value
        );

        // Update allFaceInfo from maxPointCount for all points to check
        // (usually on boundary faces)

        // Per face the new refinement data
        Map<refinementData> changedFacesInfo(pointsToCheck.size());

        forAll(pointsToCheck, i)
        {
            label pointi = pointsToCheck[i];

            // Loop over all cells using the point and check whether their
            // refinement level is much less than the maximum.

            const labelList& pCells = mesh_.pointCells(pointi);

            forAll(pCells, pCelli)
            {
                label celli = pCells[pCelli];

                refinementData& cellInfo = allCellInfo[celli];

                if
                (
                   !cellInfo.isRefined()
                 && (
                        maxPointCount[pointi]
                      > cellInfo.count() + maxFaceDiff*maxPointDiff
                    )
                )
                {
                    // Mark cell for refinement
                    cellInfo.count() = cellInfo.refinementCount();

                    // Insert faces of cell as seed faces.
                    const cell& cFaces = mesh_.cells()[celli];

                    forAll(cFaces, cFacei)
                    {
                        label facei = cFaces[cFacei];

                        refinementData faceData;
                        faceData.updateFace
                        (
                            mesh_,
                            facei,
                            celli,
                            cellInfo,
                            FaceCellWave<refinementData, int>::propagationTol(),
                            dummyTrackData
                        );

                        if (faceData.count() > allFaceInfo[facei].count())
                        {
                            changedFacesInfo.insert(facei, faceData);
                        }
                    }
                }
            }
        }

        label nChanged = changedFacesInfo.size();
        reduce(nChanged, sumOp<label>());

        if (nChanged == 0)
        {
            break;
        }


        // Transfer into seedFaces, seedFacesInfo
        seedFaces.setCapacity(changedFacesInfo.size());
        seedFacesInfo.setCapacity(changedFacesInfo.size());

        forAllConstIter(Map<refinementData>, changedFacesInfo, iter)
        {
            seedFaces.append(iter.key());
            seedFacesInfo.append(iter());
        }
    }


    if (debug)
    {
        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            label own = mesh_.faceOwner()[facei];
            label ownLevel =
                cellLevel_[own]
              + (allCellInfo[own].isRefined() ? 1 : 0);

            label nei = mesh_.faceNeighbour()[facei];
            label neiLevel =
                cellLevel_[nei]
              + (allCellInfo[nei].isRefined() ? 1 : 0);

            if (mag(ownLevel-neiLevel) > 1)
            {
                dumpCell(own);
                dumpCell(nei);
                FatalErrorInFunction
                    << "cell:" << own
                    << " current level:" << cellLevel_[own]
                    << " current refData:" << allCellInfo[own]
                    << " level after refinement:" << ownLevel
                    << nl
                    << "neighbour cell:" << nei
                    << " current level:" << cellLevel_[nei]
                    << " current refData:" << allCellInfo[nei]
                    << " level after refinement:" << neiLevel
                    << nl
                    << "which does not satisfy 2:1 constraints anymore." << nl
                    << "face:" << facei << " faceRefData:" << allFaceInfo[facei]
                    << abort(FatalError);
            }
        }


        // Coupled faces. Swap owner level to get neighbouring cell level.
        // (only boundary faces of neiLevel used)

        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
        labelList neiCount(mesh_.nFaces()-mesh_.nInternalFaces());
        labelList neiRefCount(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(neiLevel, i)
        {
            label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];
            neiLevel[i] = cellLevel_[own];
            neiCount[i] = allCellInfo[own].count();
            neiRefCount[i] = allCellInfo[own].refinementCount();
        }

        // Swap to neighbour
        syncTools::swapBoundaryFaceList(mesh_, neiLevel);
        syncTools::swapBoundaryFaceList(mesh_, neiCount);
        syncTools::swapBoundaryFaceList(mesh_, neiRefCount);

        // Now we have neighbour value see which cells need refinement
        forAll(neiLevel, i)
        {
            label facei = i+mesh_.nInternalFaces();

            label own = mesh_.faceOwner()[facei];
            label ownLevel =
                cellLevel_[own]
              + (allCellInfo[own].isRefined() ? 1 : 0);

            label nbrLevel =
                neiLevel[i]
              + ((neiCount[i] >= neiRefCount[i]) ? 1 : 0);

            if (mag(ownLevel - nbrLevel) > 1)
            {
                dumpCell(own);
                label patchi = mesh_.boundaryMesh().whichPatch(facei);

                FatalErrorInFunction
                    << "Celllevel does not satisfy 2:1 constraint."
                    << " On coupled face "
                    << facei
                    << " refData:" << allFaceInfo[facei]
                    << " on patch " << patchi << " "
                    << mesh_.boundaryMesh()[patchi].name() << nl
                    << "owner cell " << own
                    << " current level:" << cellLevel_[own]
                    << " current count:" << allCellInfo[own].count()
                    << " current refCount:"
                    << allCellInfo[own].refinementCount()
                    << " level after refinement:" << ownLevel
                    << nl
                    << "(coupled) neighbour cell"
                    << " has current level:" << neiLevel[i]
                    << " current count:" << neiCount[i]
                    << " current refCount:" << neiRefCount[i]
                    << " level after refinement:" << nbrLevel
                    << abort(FatalError);
            }
        }
    }

    // Convert back to labelList of cells to refine.

    label nRefined = 0;

    forAll(allCellInfo, celli)
    {
        if (allCellInfo[celli].isRefined())
        {
            nRefined++;
        }
    }

    // Updated list of cells to refine
    labelList newCellsToRefine(nRefined);
    nRefined = 0;

    forAll(allCellInfo, celli)
    {
        if (allCellInfo[celli].isRefined())
        {
            newCellsToRefine[nRefined++] = celli;
        }
    }

    if (debug)
    {
        Pout<< "hexRef::consistentSlowRefinement : From "
            << cellsToRefine.size() << " to " << newCellsToRefine.size()
            << " cells to refine." << endl;
    }

    return newCellsToRefine;
}


Foam::labelList Foam::hexRef::consistentSlowRefinement2
(
    const label maxFaceDiff,
    const labelList& cellsToRefine,
    const labelList& facesToCheck
) const
{
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    if (maxFaceDiff <= 0)
    {
        FatalErrorInFunction
            << "Illegal maxFaceDiff " << maxFaceDiff << nl
            << "Value should be >= 1" << exit(FatalError);
    }

    const scalar level0Size = 2*maxFaceDiff*level0EdgeLength();


    // Bit tricky. Say we want a distance of three cells between two
    // consecutive refinement levels. This is done by using FaceCellWave to
    // transport out the 'refinement shell'. Anything inside the refinement
    // shell (given by a distance) gets marked for refinement.

    // Initial information about (distance to) cellLevel on all cells
    List<refinementDistanceData> allCellInfo(mesh_.nCells());

    // Initial information about (distance to) cellLevel on all faces
    List<refinementDistanceData> allFaceInfo(mesh_.nFaces());

    // Dummy additional info for FaceCellWave
    int dummyTrackData = 0;


    // Mark cells with wanted refinement level
    forAll(cellsToRefine, i)
    {
        label celli = cellsToRefine[i];

        allCellInfo[celli] = refinementDistanceData
        (
            level0Size,
            mesh_.cellCentres()[celli],
            cellLevel_[celli]+1             // wanted refinement
        );
    }
    // Mark all others with existing refinement level
    forAll(allCellInfo, celli)
    {
        if (!allCellInfo[celli].valid(dummyTrackData))
        {
            allCellInfo[celli] = refinementDistanceData
            (
                level0Size,
                mesh_.cellCentres()[celli],
                cellLevel_[celli]           // wanted refinement
            );
        }
    }


    // Labels of seed faces
    DynamicList<label> seedFaces(mesh_.nFaces()/100);
    // refinementLevel data on seed faces
    DynamicList<refinementDistanceData> seedFacesInfo(mesh_.nFaces()/100);

    const pointField& cc = mesh_.cellCentres();

    forAll(facesToCheck, i)
    {
        label facei = facesToCheck[i];

        if (allFaceInfo[facei].valid(dummyTrackData))
        {
            // Can only occur if face has already gone through loop below.
            FatalErrorInFunction
                << "Argument facesToCheck seems to have duplicate entries!"
                << endl
                << "face:" << facei << " occurs at positions "
                << findIndices(facesToCheck, facei)
                << abort(FatalError);
        }

        label own = faceOwner[facei];

        label ownLevel =
        (
            allCellInfo[own].valid(dummyTrackData)
          ? allCellInfo[own].originLevel()
          : cellLevel_[own]
        );

        if (!mesh_.isInternalFace(facei))
        {
            // Do as if boundary face would have neighbour with one higher
            // refinement level.
            const point& fc = mesh_.faceCentres()[facei];

            refinementDistanceData neiData
            (
                level0Size,
                2*fc - cc[own],    // est'd cell centre
                ownLevel+1
            );

            allFaceInfo[facei].updateFace
            (
                mesh_,
                facei,
                own,        // not used (should be nei)
                neiData,
                FaceCellWave<refinementDistanceData, int>::propagationTol(),
                dummyTrackData
            );
        }
        else
        {
            label nei = faceNeighbour[facei];

            label neiLevel =
            (
                allCellInfo[nei].valid(dummyTrackData)
              ? allCellInfo[nei].originLevel()
              : cellLevel_[nei]
            );

            if (ownLevel == neiLevel)
            {
                // Fake as if nei>own or own>nei (whichever one 'wins')
                allFaceInfo[facei].updateFace
                (
                    mesh_,
                    facei,
                    nei,
                    refinementDistanceData(level0Size, cc[nei], neiLevel+1),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
                allFaceInfo[facei].updateFace
                (
                    mesh_,
                    facei,
                    own,
                    refinementDistanceData(level0Size, cc[own], ownLevel+1),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
            }
            else
            {
                // Difference in level anyway.
                allFaceInfo[facei].updateFace
                (
                    mesh_,
                    facei,
                    nei,
                    refinementDistanceData(level0Size, cc[nei], neiLevel),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
                allFaceInfo[facei].updateFace
                (
                    mesh_,
                    facei,
                    own,
                    refinementDistanceData(level0Size, cc[own], ownLevel),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
            }
        }
        seedFaces.append(facei);
        seedFacesInfo.append(allFaceInfo[facei]);
    }


    // Create some initial seeds to start walking from. This is only if there
    // are no facesToCheck.
    // Just seed with all faces inbetween different refinement levels for now
    forAll(faceNeighbour, facei)
    {
        // Check if face already handled in loop above
        if (!allFaceInfo[facei].valid(dummyTrackData))
        {
            label own = faceOwner[facei];

            label ownLevel =
            (
                allCellInfo[own].valid(dummyTrackData)
              ? allCellInfo[own].originLevel()
              : cellLevel_[own]
            );

            label nei = faceNeighbour[facei];

            label neiLevel =
            (
                allCellInfo[nei].valid(dummyTrackData)
              ? allCellInfo[nei].originLevel()
              : cellLevel_[nei]
            );

            if (ownLevel > neiLevel)
            {
                // Set face to owner data. (since face not yet would be copy)
                seedFaces.append(facei);
                allFaceInfo[facei].updateFace
                (
                    mesh_,
                    facei,
                    own,
                    refinementDistanceData(level0Size, cc[own], ownLevel),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
                seedFacesInfo.append(allFaceInfo[facei]);
            }
            else if (neiLevel > ownLevel)
            {
                seedFaces.append(facei);
                allFaceInfo[facei].updateFace
                (
                    mesh_,
                    facei,
                    nei,
                    refinementDistanceData(level0Size, cc[nei], neiLevel),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
                seedFacesInfo.append(allFaceInfo[facei]);
            }
        }
    }

    seedFaces.shrink();
    seedFacesInfo.shrink();

    // face-cell-face transport engine
    FaceCellWave<refinementDistanceData, int> levelCalc
    (
        mesh_,
        seedFaces,
        seedFacesInfo,
        allFaceInfo,
        allCellInfo,
        mesh_.globalData().nTotalCells()+1,
        dummyTrackData
    );


    //if (debug)
    //{
    //    // Dump wanted level
    //    volScalarField wantedLevel
    //    (
    //        IOobject
    //        (
    //            "wantedLevel",
    //            fMesh.time().timeName(),
    //            fMesh,
    //            IOobject::NO_READ,
    //            IOobject::NO_WRITE,
    //            false
    //        ),
    //        fMesh,
    //        dimensionedScalar("zero", dimless, 0)
    //    );
    //
    //    forAll(wantedLevel, celli)
    //    {
    //        wantedLevel[celli] = allCellInfo[celli].wantedLevel(cc[celli]);
    //    }
    //
    //    Pout<< "Writing " << wantedLevel.objectPath() << endl;
    //    wantedLevel.write();
    //}


    // Convert back to labelList of cells to refine.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // 1. Force original refinement cells to be picked up by setting the
    // originLevel of input cells to be a very large level (but within range
    // of 1<< shift inside refinementDistanceData::wantedLevel)
    forAll(cellsToRefine, i)
    {
        label celli = cellsToRefine[i];

        allCellInfo[celli].originLevel() = sizeof(label)*8-2;
        allCellInfo[celli].origin() = cc[celli];
    }

    // 2. Extend to 2:1. I don't understand yet why this is not done
    // 2. Extend to 2:1. For non-cube cells the scalar distance does not work
    // so make sure it at least provides 2:1.
    PackedBoolList refineCell(mesh_.nCells());
    forAll(allCellInfo, celli)
    {
        label wanted = allCellInfo[celli].wantedLevel(cc[celli]);

        if (wanted > cellLevel_[celli]+1)
        {
            refineCell.set(celli);
        }
    }
    faceConsistentRefinement(true, refineCell);

    while (true)
    {
        label nChanged = faceConsistentRefinement(true, refineCell);

        reduce(nChanged, sumOp<label>());

        if (debug)
        {
            Pout<< "hexRef::consistentSlowRefinement2 : Changed " << nChanged
                << " refinement levels due to 2:1 conflicts."
                << endl;
        }

        if (nChanged == 0)
        {
            break;
        }
    }

    // 3. Convert back to labelList.
    label nRefined = 0;

    forAll(refineCell, celli)
    {
//        if (refineCell.get(celli))
        if (refineCell[celli])
        {
            nRefined++;
        }
    }

    labelList newCellsToRefine(nRefined);
    nRefined = 0;

    forAll(refineCell, celli)
    {
//        if (refineCell.get(celli))
        if (refineCell[celli])
        {
            newCellsToRefine[nRefined++] = celli;
        }
    }

    if (debug)
    {
        Pout<< "hexRef::consistentSlowRefinement2 : From "
            << cellsToRefine.size() << " to " << newCellsToRefine.size()
            << " cells to refine." << endl;

        // Check that newCellsToRefine obeys at least 2:1.

        {
            cellSet cellsIn(mesh_, "cellsToRefineIn", cellsToRefine);
            Pout<< "hexRef::consistentSlowRefinement2 : writing "
                << cellsIn.size() << " to cellSet "
                << cellsIn.objectPath() << endl;
            cellsIn.write();
        }
        {
            cellSet cellsOut(mesh_, "cellsToRefineOut", newCellsToRefine);
            Pout<< "hexRef::consistentSlowRefinement2 : writing "
                << cellsOut.size() << " to cellSet "
                << cellsOut.objectPath() << endl;
            cellsOut.write();
        }

        // Extend to 2:1
        PackedBoolList refineCell(mesh_.nCells());
        forAll(newCellsToRefine, i)
        {
            refineCell.set(newCellsToRefine[i]);
        }
        const PackedBoolList savedRefineCell(refineCell);

        label nChanged = faceConsistentRefinement(true, refineCell);

        {
            cellSet cellsOut2
            (
                mesh_, "cellsToRefineOut2", newCellsToRefine.size()
            );
            forAll(refineCell, celli)
            {
                if (refineCell.get(celli))
                {
                    cellsOut2.insert(celli);
                }
            }
            Pout<< "hexRef::consistentSlowRefinement2 : writing "
                << cellsOut2.size() << " to cellSet "
                << cellsOut2.objectPath() << endl;
            cellsOut2.write();
        }

        if (nChanged > 0)
        {
            forAll(refineCell, celli)
            {
                if (refineCell.get(celli) && !savedRefineCell.get(celli))
                {
                    dumpCell(celli);
                    FatalErrorInFunction
                        << "Cell:" << celli << " cc:"
                        << mesh_.cellCentres()[celli]
                        << " was not marked for refinement but does not obey"
                        << " 2:1 constraints."
                        << abort(FatalError);
                }
            }
        }
    }

    return newCellsToRefine;
}


void Foam::hexRef::storeData
(
    const labelList& pointsToStore,
    const labelList& facesToStore,
    const labelList& cellsToStore
)
{
    savedPointLevel_.resize(2*pointsToStore.size());
    forAll(pointsToStore, i)
    {
        label pointi = pointsToStore[i];
        savedPointLevel_.insert(pointi, pointLevel_[pointi]);
    }

    savedCellLevel_.resize(2*cellsToStore.size());
    forAll(cellsToStore, i)
    {
        label celli = cellsToStore[i];
        savedCellLevel_.insert(celli, cellLevel_[celli]);
    }
}


// Gets called after the mesh change. setRefinement will already have made
// sure the pointLevel_ and cellLevel_ are the size of the new mesh so we
// only need to account for reordering.
void Foam::hexRef::updateMesh(const mapPolyMesh& map)
{
    Map<label> dummyMap(0);

    updateMesh(map, dummyMap, dummyMap, dummyMap);
}


// Gets called after the mesh change. setRefinement will already have made
// sure the pointLevel_ and cellLevel_ are the size of the new mesh so we
// only need to account for reordering.
void Foam::hexRef::updateMesh
(
    const mapPolyMesh& map,
    const Map<label>& pointsToRestore,
    const Map<label>& facesToRestore,
    const Map<label>& cellsToRestore
)
{
    // Update celllevel
    if (debug)
    {
        Pout<< "hexRef::updateMesh :"
            << " Updating various lists"
            << endl;
    }

    {
        const labelList& reverseCellMap = map.reverseCellMap();

        if (debug)
        {
            Pout<< "hexRef::updateMesh :"
                << " reverseCellMap:" << map.reverseCellMap().size()
                << " cellMap:" << map.cellMap().size()
                << " nCells:" << mesh_.nCells()
                << " nOldCells:" << map.nOldCells()
                << " cellLevel_:" << cellLevel_.size()
                << " reversePointMap:" << map.reversePointMap().size()
                << " pointMap:" << map.pointMap().size()
                << " nPoints:" << mesh_.nPoints()
                << " nOldPoints:" << map.nOldPoints()
                << " pointLevel_:" << pointLevel_.size()
                << endl;
        }

        if (reverseCellMap.size() == cellLevel_.size())
        {
            // Assume it is after hexRef that this routine is called.
            // Just account for reordering. We cannot use cellMap since
            // then cells created from cells would get cellLevel_ of
            // cell they were created from.
            reorder(reverseCellMap, mesh_.nCells(), -1, cellLevel_);
        }
        else
        {
            // Map data
            const labelList& cellMap = map.cellMap();

            labelList newCellLevel(cellMap.size());
            forAll(cellMap, newCelli)
            {
                label oldCelli = cellMap[newCelli];

                if (oldCelli == -1)
                {
                    newCellLevel[newCelli] = -1;
                }
                else
                {
                    newCellLevel[newCelli] = cellLevel_[oldCelli];
                }
            }
            cellLevel_.transfer(newCellLevel);
        }

        // See if any cells to restore. This will be for some new cells
        // the corresponding old cell.
        forAllConstIter(Map<label>, cellsToRestore, iter)
        {
            label newCelli = iter.key();
            label storedCelli = iter();

            Map<label>::iterator fnd = savedCellLevel_.find(storedCelli);

            if (fnd == savedCellLevel_.end())
            {
                FatalErrorInFunction
                    << "Problem : trying to restore old value for new cell "
                    << newCelli << " but cannot find old cell " << storedCelli
                    << " in map of stored values " << savedCellLevel_
                    << abort(FatalError);
            }
            cellLevel_[newCelli] = fnd();
        }

        //if (findIndex(cellLevel_, -1) != -1)
        //{
        //    WarningInFunction
        //        << "Problem : "
        //        << "cellLevel_ contains illegal value -1 after mapping
        //        << " at cell " << findIndex(cellLevel_, -1) << endl
        //        << "This means that another program has inflated cells"
        //        << " (created cells out-of-nothing) and hence we don't know"
        //        << " their cell level. Continuing with illegal value."
        //        << abort(FatalError);
        //}
    }


    // Update pointlevel
    {
        const labelList& reversePointMap = map.reversePointMap();

        if (reversePointMap.size() == pointLevel_.size())
        {
            // Assume it is after hexRef that this routine is called.
            reorder(reversePointMap, mesh_.nPoints(), -1,  pointLevel_);
        }
        else
        {
            // Map data
            const labelList& pointMap = map.pointMap();

            labelList newPointLevel(pointMap.size());

            forAll(pointMap, newPointi)
            {
                label oldPointi = pointMap[newPointi];

                if (oldPointi == -1)
                {
                    //FatalErrorInFunction
                    //    << "Problem : point " << newPointi
                    //    << " at " << mesh_.points()[newPointi]
                    //    << " does not originate from another point"
                    //    << " (i.e. is inflated)." << nl
                    //    << "Hence we cannot determine the new pointLevel"
                    //    << " for it." << abort(FatalError);
                    newPointLevel[newPointi] = -1;
                }
                else
                {
                    newPointLevel[newPointi] = pointLevel_[oldPointi];
                }
            }
            pointLevel_.transfer(newPointLevel);
        }

        // See if any points to restore. This will be for some new points
        // the corresponding old point (the one from the call to storeData)
        forAllConstIter(Map<label>, pointsToRestore, iter)
        {
            label newPointi = iter.key();
            label storedPointi = iter();

            Map<label>::iterator fnd = savedPointLevel_.find(storedPointi);

            if (fnd == savedPointLevel_.end())
            {
                FatalErrorInFunction
                    << "Problem : trying to restore old value for new point "
                    << newPointi << " but cannot find old point "
                    << storedPointi
                    << " in map of stored values " << savedPointLevel_
                    << abort(FatalError);
            }
            pointLevel_[newPointi] = fnd();
        }

        //if (findIndex(pointLevel_, -1) != -1)
        //{
        //    WarningInFunction
        //        << "Problem : "
        //        << "pointLevel_ contains illegal value -1 after mapping"
        //        << " at point" << findIndex(pointLevel_, -1) << endl
        //        << "This means that another program has inflated points"
        //        << " (created points out-of-nothing) and hence we don't know"
        //        << " their point level. Continuing with illegal value."
        //        //<< abort(FatalError);
        //}
    }

    // Update refinement tree
    if (history_.active())
    {
        history_.updateMesh(map);
    }

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    // Update face removal engine
    faceRemover_.updateMesh(map);

    // Clear cell shapes
    cellShapesPtr_.clear();
}


// Gets called after mesh subsetting. Maps are from new back to old.
void Foam::hexRef::subset
(
    const labelList& pointMap,
    const labelList& faceMap,
    const labelList& cellMap
)
{
    // Update celllevel
    if (debug)
    {
        Pout<< "hexRef::subset :"
            << " Updating various lists"
            << endl;
    }

    if (history_.active())
    {
        WarningInFunction
            << "Subsetting will not work in combination with unrefinement."
            << nl
            << "Proceed at your own risk." << endl;
    }


    // Update celllevel
    {
        labelList newCellLevel(cellMap.size());

        forAll(cellMap, newCelli)
        {
            newCellLevel[newCelli] = cellLevel_[cellMap[newCelli]];
        }

        cellLevel_.transfer(newCellLevel);

        if (findIndex(cellLevel_, -1) != -1)
        {
            FatalErrorInFunction
                << "Problem : "
                << "cellLevel_ contains illegal value -1 after mapping:"
                << cellLevel_
                << abort(FatalError);
        }
    }

    // Update pointlevel
    {
        labelList newPointLevel(pointMap.size());

        forAll(pointMap, newPointi)
        {
            newPointLevel[newPointi] = pointLevel_[pointMap[newPointi]];
        }

        pointLevel_.transfer(newPointLevel);

        if (findIndex(pointLevel_, -1) != -1)
        {
            FatalErrorInFunction
                << "Problem : "
                << "pointLevel_ contains illegal value -1 after mapping:"
                << pointLevel_
                << abort(FatalError);
        }
    }

    // Update refinement tree
    if (history_.active())
    {
        history_.subset(pointMap, faceMap, cellMap);
    }

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    // Nothing needs doing to faceRemover.
    //faceRemover_.subset(pointMap, faceMap, cellMap);

    // Clear cell shapes
    cellShapesPtr_.clear();
}


// Gets called after the mesh distribution
void Foam::hexRef::distribute(const mapDistributePolyMesh& map)
{
    if (debug)
    {
        Pout<< "hexRef::distribute :"
            << " Updating various lists"
            << endl;
    }

    // Update celllevel
    map.distributeCellData(cellLevel_);
    // Update pointlevel
    map.distributePointData(pointLevel_);

    // Update refinement tree
    if (history_.active())
    {
        history_.distribute(map);
    }

    // Update face removal engine
    faceRemover_.distribute(map);

    // Clear cell shapes
    cellShapesPtr_.clear();
}


void Foam::hexRef::checkMesh() const
{
    const scalar smallDim = 1e-6 * mesh_.bounds().mag();

    if (debug)
    {
        Pout<< "hexRef::checkMesh : Using matching tolerance smallDim:"
            << smallDim << endl;
    }

    // Check owner on coupled faces.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // There should be only one coupled face between two cells. Why? Since
    // otherwise mesh redistribution might cause multiple faces between two
    // cells
    {
        labelList nei(mesh_.nFaces()-mesh_.nInternalFaces());
        forAll(nei, i)
        {
            nei[i] = mesh_.faceOwner()[i+mesh_.nInternalFaces()];
        }

        // Replace data on coupled patches with their neighbour ones.
        syncTools::swapBoundaryFaceList(mesh_, nei);

        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];

            if (pp.coupled())
            {
                // Check how many faces between owner and neighbour. Should
                // be only one.
                HashTable<label, labelPair, labelPair::Hash<>>
                    cellToFace(2*pp.size());

                label facei = pp.start();

                forAll(pp, i)
                {
                    label own = mesh_.faceOwner()[facei];
                    label bFacei = facei-mesh_.nInternalFaces();

                    if (!cellToFace.insert(labelPair(own, nei[bFacei]), facei))
                    {
                        dumpCell(own);
                        FatalErrorInFunction
                            << "Faces do not seem to be correct across coupled"
                            << " boundaries" << endl
                            << "Coupled face " << facei
                            << " between owner " << own
                            << " on patch " << pp.name()
                            << " and coupled neighbour " << nei[bFacei]
                            << " has two faces connected to it:"
                            << facei << " and "
                            << cellToFace[labelPair(own, nei[bFacei])]
                            << abort(FatalError);
                    }

                    facei++;
                }
            }
        }
    }

    // Check face areas.
    // ~~~~~~~~~~~~~~~~~

    {
        scalarField neiFaceAreas(mesh_.nFaces()-mesh_.nInternalFaces());
        forAll(neiFaceAreas, i)
        {
            neiFaceAreas[i] = mag(mesh_.faceAreas()[i+mesh_.nInternalFaces()]);
        }

        // Replace data on coupled patches with their neighbour ones.
        syncTools::swapBoundaryFaceList(mesh_, neiFaceAreas);

        forAll(neiFaceAreas, i)
        {
            label facei = i+mesh_.nInternalFaces();

            const scalar magArea = mag(mesh_.faceAreas()[facei]);

            if (mag(magArea - neiFaceAreas[i]) > smallDim)
            {
                const face& f = mesh_.faces()[facei];
                label patchi = mesh_.boundaryMesh().whichPatch(facei);

                dumpCell(mesh_.faceOwner()[facei]);

                FatalErrorInFunction
                    << "Faces do not seem to be correct across coupled"
                    << " boundaries" << endl
                    << "Coupled face " << facei
                    << " on patch " << patchi
                    << " " << mesh_.boundaryMesh()[patchi].name()
                    << " coords:" << UIndirectList<point>(mesh_.points(), f)()
                    << " has face area:" << magArea
                    << " (coupled) neighbour face area differs:"
                    << neiFaceAreas[i]
                    << " to within tolerance " << smallDim
                    << abort(FatalError);
            }
        }
    }


    // Check number of points on faces.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    {
        labelList nVerts(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(nVerts, i)
        {
            nVerts[i] = mesh_.faces()[i+mesh_.nInternalFaces()].size();
        }

        // Replace data on coupled patches with their neighbour ones.
        syncTools::swapBoundaryFaceList(mesh_, nVerts);

        forAll(nVerts, i)
        {
            label facei = i+mesh_.nInternalFaces();

            const face& f = mesh_.faces()[facei];

            if (f.size() != nVerts[i])
            {
                dumpCell(mesh_.faceOwner()[facei]);

                label patchi = mesh_.boundaryMesh().whichPatch(facei);

                FatalErrorInFunction
                    << "Faces do not seem to be correct across coupled"
                    << " boundaries" << endl
                    << "Coupled face " << facei
                    << " on patch " << patchi
                    << " " << mesh_.boundaryMesh()[patchi].name()
                    << " coords:" << UIndirectList<point>(mesh_.points(), f)()
                    << " has size:" << f.size()
                    << " (coupled) neighbour face has size:"
                    << nVerts[i]
                    << abort(FatalError);
            }
        }
    }


    // Check points of face
    // ~~~~~~~~~~~~~~~~~~~~
    {
        // Anchor points.
        pointField anchorPoints(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(anchorPoints, i)
        {
            label facei = i+mesh_.nInternalFaces();
            const point& fc = mesh_.faceCentres()[facei];
            const face& f = mesh_.faces()[facei];
            const vector anchorVec(mesh_.points()[f[0]] - fc);

            anchorPoints[i] = anchorVec;
        }

        // Replace data on coupled patches with their neighbour ones. Apply
        // rotation transformation (but not separation since is relative vector
        // to point on same face.
        syncTools::swapBoundaryFaceList(mesh_, anchorPoints);

        forAll(anchorPoints, i)
        {
            label facei = i+mesh_.nInternalFaces();
            const point& fc = mesh_.faceCentres()[facei];
            const face& f = mesh_.faces()[facei];
            const vector anchorVec(mesh_.points()[f[0]] - fc);

            if (mag(anchorVec - anchorPoints[i]) > smallDim)
            {
                dumpCell(mesh_.faceOwner()[facei]);

                label patchi = mesh_.boundaryMesh().whichPatch(facei);

                FatalErrorInFunction
                    << "Faces do not seem to be correct across coupled"
                    << " boundaries" << endl
                    << "Coupled face " << facei
                    << " on patch " << patchi
                    << " " << mesh_.boundaryMesh()[patchi].name()
                    << " coords:" << UIndirectList<point>(mesh_.points(), f)()
                    << " has anchor vector:" << anchorVec
                    << " (coupled) neighbour face anchor vector differs:"
                    << anchorPoints[i]
                    << " to within tolerance " << smallDim
                    << abort(FatalError);
            }
        }
    }

    if (debug)
    {
        Pout<< "hexRef::checkMesh : Returning" << endl;
    }
}


void Foam::hexRef::checkRefinementLevels
(
    const label maxPointDiff,
    const labelList& pointsToCheck
) const
{
    if (debug)
    {
        Pout<< "hexRef::checkRefinementLevels :"
            << " Checking 2:1 refinement level" << endl;
    }

    if
    (
        cellLevel_.size() != mesh_.nCells()
     || pointLevel_.size() != mesh_.nPoints()
    )
    {
        FatalErrorInFunction
            << "cellLevel size should be number of cells"
            << " and pointLevel size should be number of points."<< nl
            << "cellLevel:" << cellLevel_.size()
            << " mesh.nCells():" << mesh_.nCells() << nl
            << "pointLevel:" << pointLevel_.size()
            << " mesh.nPoints():" << mesh_.nPoints()
            << abort(FatalError);
    }


    // Check 2:1 consistency.
    // ~~~~~~~~~~~~~~~~~~~~~~

    {
        // Internal faces.
        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            label own = mesh_.faceOwner()[facei];
            label nei = mesh_.faceNeighbour()[facei];

            if (mag(cellLevel_[own] - cellLevel_[nei]) > 1)
            {
                dumpCell(own);
                dumpCell(nei);

                FatalErrorInFunction
                    << "Celllevel does not satisfy 2:1 constraint." << nl
                    << "On face " << facei << " owner cell " << own
                    << " has refinement " << cellLevel_[own]
                    << " neighbour cell " << nei << " has refinement "
                    << cellLevel_[nei]
                    << abort(FatalError);
            }
        }

        // Coupled faces. Get neighbouring value
        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(neiLevel, i)
        {
            label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];

            neiLevel[i] = cellLevel_[own];
        }

        // No separation
        syncTools::swapBoundaryFaceList(mesh_, neiLevel);

        forAll(neiLevel, i)
        {
            label facei = i+mesh_.nInternalFaces();

            label own = mesh_.faceOwner()[facei];

            if (mag(cellLevel_[own] - neiLevel[i]) > 1)
            {
                dumpCell(own);

                label patchi = mesh_.boundaryMesh().whichPatch(facei);

                FatalErrorInFunction
                    << "Celllevel does not satisfy 2:1 constraint."
                    << " On coupled face " << facei
                    << " on patch " << patchi << " "
                    << mesh_.boundaryMesh()[patchi].name()
                    << " owner cell " << own << " has refinement "
                    << cellLevel_[own]
                    << " (coupled) neighbour cell has refinement "
                    << neiLevel[i]
                    << abort(FatalError);
            }
        }
    }


    // Check pointLevel is synchronized
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    {
        labelList syncPointLevel(pointLevel_);

        // Get min level
        syncTools::syncPointList
        (
            mesh_,
            syncPointLevel,
            minEqOp<label>(),
            labelMax
        );


        forAll(syncPointLevel, pointi)
        {
            if (pointLevel_[pointi] != syncPointLevel[pointi])
            {
                FatalErrorInFunction
                    << "PointLevel is not consistent across coupled patches."
                    << endl
                    << "point:" << pointi << " coord:" << mesh_.points()[pointi]
                    << " has level " << pointLevel_[pointi]
                    << " whereas the coupled point has level "
                    << syncPointLevel[pointi]
                    << abort(FatalError);
            }
        }
    }


    // Check 2:1 across points (instead of faces)
    if (maxPointDiff != -1)
    {
        // Determine per point the max cell level.
        labelList maxPointLevel(mesh_.nPoints(), 0);

        forAll(maxPointLevel, pointi)
        {
            const labelList& pCells = mesh_.pointCells(pointi);

            label& pLevel = maxPointLevel[pointi];

            forAll(pCells, i)
            {
                pLevel = max(pLevel, cellLevel_[pCells[i]]);
            }
        }

        // Sync maxPointLevel to neighbour
        syncTools::syncPointList
        (
            mesh_,
            maxPointLevel,
            maxEqOp<label>(),
            labelMin            // null value
        );

        // Check 2:1 across boundary points
        forAll(pointsToCheck, i)
        {
            label pointi = pointsToCheck[i];

            const labelList& pCells = mesh_.pointCells(pointi);

            forAll(pCells, i)
            {
                label celli = pCells[i];

                if
                (
                    mag(cellLevel_[celli]-maxPointLevel[pointi])
                  > maxPointDiff
                )
                {
                    dumpCell(celli);

                    FatalErrorInFunction
                        << "Too big a difference between"
                        << " point-connected cells." << nl
                        << "cell:" << celli
                        << " cellLevel:" << cellLevel_[celli]
                        << " uses point:" << pointi
                        << " coord:" << mesh_.points()[pointi]
                        << " which is also used by a cell with level:"
                        << maxPointLevel[pointi]
                        << abort(FatalError);
                }
            }
        }
    }


    //- Gives problems after first splitting off inside mesher.
    //// Hanging points
    //{
    //    // Any patches with points having only two edges.
    //
    //    boolList isHangingPoint(mesh_.nPoints(), false);
    //
    //    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    //
    //    forAll(patches, patchi)
    //    {
    //        const polyPatch& pp = patches[patchi];
    //
    //        const labelList& meshPoints = pp.meshPoints();
    //
    //        forAll(meshPoints, i)
    //        {
    //            label pointi = meshPoints[i];
    //
    //            const labelList& pEdges = mesh_.pointEdges()[pointi];
    //
    //            if (pEdges.size() == 2)
    //            {
    //                isHangingPoint[pointi] = true;
    //            }
    //        }
    //    }
    //
    //    syncTools::syncPointList
    //    (
    //        mesh_,
    //        isHangingPoint,
    //        andEqOp<bool>(),        // only if all decide it is hanging point
    //        true,                   // null
    //        false                   // no separation
    //    );
    //
    //    //OFstream str(mesh_.time().path()/"hangingPoints.obj");
    //
    //    label nHanging = 0;
    //
    //    forAll(isHangingPoint, pointi)
    //    {
    //        if (isHangingPoint[pointi])
    //        {
    //            nHanging++;
    //
    //            Pout<< "Hanging boundary point " << pointi
    //                << " at " << mesh_.points()[pointi]
    //                << endl;
    //            //meshTools::writeOBJ(str, mesh_.points()[pointi]);
    //        }
    //    }
    //
    //    if (returnReduce(nHanging, sumOp<label>()) > 0)
    //    {
    //        FatalErrorInFunction
    //            << "Detected a point used by two edges only (hanging point)"
    //            << nl << "This is not allowed"
    //            << abort(FatalError);
    //    }
    //}
}


const Foam::cellShapeList& Foam::hexRef::cellShapes() const
{
    if (cellShapesPtr_.empty())
    {
        if (debug)
        {
            Pout<< "hexRef::cellShapes() : calculating splitHex cellShapes."
                << " cellLevel:" << cellLevel_.size()
                << " pointLevel:" << pointLevel_.size()
                << endl;
        }

        const cellShapeList& meshShapes = mesh_.cellShapes();
        cellShapesPtr_.reset(new cellShapeList(meshShapes));

        label nSplitHex = 0;
        label nUnrecognised = 0;

        forAll(cellLevel_, celli)
        {
            if (meshShapes[celli].model().index() == 0)
            {
                label level = cellLevel_[celli];

                // Return true if we've found 6 quads
                DynamicList<face> quads;
                bool haveQuads = matchHexShape
                (
                    celli,
                    level,
                    quads
                );

                if (haveQuads)
                {
                    faceList faces(move(quads));
                    cellShapesPtr_()[celli] = degenerateMatcher::match(faces);
                    nSplitHex++;
                }
                else
                {
                    nUnrecognised++;
                }
            }
        }
        if (debug)
        {
            Pout<< "hexRef::cellShapes() :"
                << " nCells:" << mesh_.nCells() << " of which" << nl
                << "    primitive:" << (mesh_.nCells()-nSplitHex-nUnrecognised)
                << nl
                << "    split-hex:" << nSplitHex << nl
                << "    poly     :" << nUnrecognised << nl
                << endl;
        }
    }
    return cellShapesPtr_();
}



// Write refinement to polyMesh directory.
bool Foam::hexRef::write() const
{
    bool writeOk =
        cellLevel_.write()
     && pointLevel_.write()
     && level0Edge_.write();

    if (history_.active())
    {
        writeOk = writeOk && history_.write();
    }

    return writeOk;
}




// ************************************************************************* //
