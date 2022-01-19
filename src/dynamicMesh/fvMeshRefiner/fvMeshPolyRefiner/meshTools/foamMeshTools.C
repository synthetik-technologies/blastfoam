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

#include "foamMeshTools.H"
#include "polyMesh.H"
#include "hexMatcher.H"
#include "faceZone.H"
#include "syncTools.H"

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
/*

// Does anyone have couples? Since meshes might have 0 cells and 0 proc
// boundaries need to reduce this info.
bool Foam::foamSyncTools::hasCouples(const polyBoundaryMesh& patches)
{
    bool hasAnyCouples = false;

    forAll(patches, patchI)
    {
        if (patches[patchI].coupled())
        {
            hasAnyCouples = true;
            break;
        }
    }
    return returnReduce(hasAnyCouples, orOp<bool>());
}


void Foam::foamSyncTools::checkTransform
(
    const coupledPolyPatch& pp,
    const bool applySeparation
)
{
    if (!pp.parallel() && pp.forwardT().size() > 1)
    {
        FatalErrorInFunction
            << "Non-uniform transformation not supported for point or edge"
            << " fields." << endl
            << "Patch:" << pp.name()
            << abort(FatalError);
    }
    if (applySeparation && pp.separated() && pp.separation().size() > 1)
    {
        FatalErrorInFunction
            << "Non-uniform separation vector not supported for point or edge"
            << " fields." << endl
            << "Patch:" << pp.name()
            << abort(FatalError);
    }
}


// Determines for every point whether it is coupled and if so sets only one.
Foam::PackedBoolList Foam::foamSyncTools::getMasterPoints(const polyMesh& mesh)
{
    PackedBoolList isMasterPoint(mesh.nPoints(), 0);
    PackedBoolList donePoint(mesh.nPoints(), 0);


    // Do multiple shared points. Min. proc is master
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const labelList& sharedPointAddr =
        mesh.globalData().sharedPointAddr();

    labelList minProc(mesh.globalData().nGlobalPoints(), labelMax);

    UIndirectList<label>(minProc, sharedPointAddr) = Pstream::myProcNo();

    Pstream::listCombineGather(minProc, minEqOp<label>());
    Pstream::listCombineScatter(minProc);

    const labelList& sharedPointLabels =
        mesh.globalData().sharedPointLabels();

    forAll(sharedPointAddr, i)
    {
        if (minProc[sharedPointAddr[i]] == Pstream::myProcNo())
        {
            isMasterPoint.set(sharedPointLabels[i], 1u);
        }
        donePoint.set(sharedPointLabels[i], 1u);
    }


    // Do other points on coupled patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchI)
    {
        if (patches[patchI].coupled())
        {
            if
            (
                Pstream::parRun()
             && isA<processorPolyPatch>(patches[patchI])
            )
            {
                const processorPolyPatch& pp =
                    refCast<const processorPolyPatch>(patches[patchI]);

                const labelList& meshPoints = pp.meshPoints();

                forAll(meshPoints, i)
                {
                    label pointI = meshPoints[i];

                    if (donePoint.get(pointI) == 0u)
                    {
                        donePoint.set(pointI, 1u);

                        if (pp.master())
                        {
                            isMasterPoint.set(pointI, 1u);
                        }
                    }
                }
            }
            else if (isA<cyclicPolyPatch>(patches[patchI]))
            {
                const cyclicPolyPatch& pp =
                    refCast<const cyclicPolyPatch>(patches[patchI]);

                const edgeList& coupledPoints = pp.coupledPoints();
                const labelList& meshPoints = pp.meshPoints();

                forAll(coupledPoints, i)
                {
                    // First one of couple points is master

                    const edge& pointPair = coupledPoints[i];
                    label p0 = meshPoints[pointPair[0]];
                    label p1 = meshPoints[pointPair[1]];

                    if (donePoint.get(p0) == 0u)
                    {
                        donePoint.set(p0, 1u);
                        isMasterPoint.set(p0, 1u);
                        donePoint.set(p1, 1u);
                    }
                }
            }
            else
            {
                FatalErrorInFunction
                    << "Cannot handle coupled patch " << patches[patchI].name()
                    << " of type " <<  patches[patchI].type()
                    << abort(FatalError);
            }
        }
    }


    // Do all other points
    // ~~~~~~~~~~~~~~~~~~~

    forAll(donePoint, pointI)
    {
        if (donePoint.get(pointI) == 0u)
        {
            donePoint.set(pointI, 1u);
            isMasterPoint.set(pointI, 1u);
        }
    }

    return isMasterPoint;
}


// Determines for every edge whether it is coupled and if so sets only one.
Foam::PackedBoolList Foam::foamSyncTools::getMasterEdges(const polyMesh& mesh)
{
    PackedBoolList isMasterEdge(mesh.nEdges(), 0);
    PackedBoolList doneEdge(mesh.nEdges(), 0);


    // Do multiple shared edges. Min. proc is master
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const labelList& sharedEdgeAddr =
        mesh.globalData().sharedEdgeAddr();

    labelList minProc(mesh.globalData().nGlobalEdges(), labelMax);

    UIndirectList<label>(minProc, sharedEdgeAddr) = Pstream::myProcNo();

    Pstream::listCombineGather(minProc, minEqOp<label>());
    Pstream::listCombineScatter(minProc);

    const labelList& sharedEdgeLabels =
        mesh.globalData().sharedEdgeLabels();

    forAll(sharedEdgeAddr, i)
    {
        if (minProc[sharedEdgeAddr[i]] == Pstream::myProcNo())
        {
            isMasterEdge.set(sharedEdgeLabels[i], 1u);
        }
        doneEdge.set(sharedEdgeLabels[i], 1u);
    }


    // Do other edges on coupled patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchI)
    {
        if (patches[patchI].coupled())
        {
            if
            (
                Pstream::parRun()
             && isA<processorPolyPatch>(patches[patchI])
            )
            {
                const processorPolyPatch& pp =
                    refCast<const processorPolyPatch>(patches[patchI]);

                const labelList& meshEdges = pp.meshEdges();

                forAll(meshEdges, i)
                {
                    label edgeI = meshEdges[i];

                    if (doneEdge.get(edgeI) == 0u)
                    {
                        doneEdge.set(edgeI, 1u);

                        if (pp.master())
                        {
                            isMasterEdge.set(edgeI, 1u);
                        }
                    }
                }
            }
            else if (isA<cyclicPolyPatch>(patches[patchI]))
            {
                const cyclicPolyPatch& pp =
                    refCast<const cyclicPolyPatch>(patches[patchI]);

                const edgeList& coupledEdges = pp.coupledEdges();
                const labelList& meshEdges = pp.meshEdges();

                forAll(coupledEdges, i)
                {
                    // First one of couple edges is master

                    const edge& edgePair = coupledEdges[i];
                    label e0 = meshEdges[edgePair[0]];
                    label e1 = meshEdges[edgePair[1]];

                    if (doneEdge.get(e0) == 0u)
                    {
                        doneEdge.set(e0, 1u);
                        isMasterEdge.set(e0, 1u);
                        doneEdge.set(e1, 1u);
                    }
                }
            }
            else
            {
                FatalErrorInFunction
                    << "Cannot handle coupled patch " << patches[patchI].name()
                    << " of type " <<  patches[patchI].type()
                    << abort(FatalError);
            }
        }
    }


    // Do all other edges
    // ~~~~~~~~~~~~~~~~~~

    forAll(doneEdge, edgeI)
    {
        if (doneEdge.get(edgeI) == 0u)
        {
            doneEdge.set(edgeI, 1u);
            isMasterEdge.set(edgeI, 1u);
        }
    }

    return isMasterEdge;
}


// Determines for every face whether it is coupled and if so sets only one.
Foam::PackedBoolList Foam::foamSyncTools::getMasterFaces(const polyMesh& mesh)
{
    PackedBoolList isMasterFace(mesh.nFaces(), 1);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchI)
    {
        if (patches[patchI].coupled())
        {
            if (Pstream::parRun() && isA<processorPolyPatch>(patches[patchI]))
            {
                const processorPolyPatch& pp =
                    refCast<const processorPolyPatch>(patches[patchI]);

                if (!pp.master())
                {
                    forAll(pp, i)
                    {
                        isMasterFace.set(pp.start()+i, 0);
                    }
                }
            }
            else if (isA<cyclicPolyPatch>(patches[patchI]))
            {
                const cyclicPolyPatch& pp =
                    refCast<const cyclicPolyPatch>(patches[patchI]);

                for (label i = pp.size()/2; i < pp.size(); i++)
                {
                    isMasterFace.set(pp.start()+i, 0);
                }
            }
            else
            {
                FatalErrorInFunction
                    << "Cannot handle coupled patch " << patches[patchI].name()
                    << " of type " <<  patches[patchI].type()
                    << abort(FatalError);
            }
        }
    }

    return isMasterFace;
}


template <>
void Foam::syncTools::foamSyncTools
(
    const vectorField& separation,
    UList<vector>& field
)
{
    if (separation.size() == 1)
    {
        // Single value for all.

        forAll(field, i)
        {
            field[i] += separation[0];
        }
    }
    else if (separation.size() == field.size())
    {
        forAll(field, i)
        {
            field[i] += separation[i];
        }
    }
    else
    {
        FatalErrorInFunction
            << "Sizes of field and transformation not equal. field:"
            << field.size() << " transformation:" << separation.size()
            << abort(FatalError);
    }
}


template <>
void Foam::foamSyncTools::separateList
(
    const vectorField& separation,
    Map<vector>& field
)
{
    if (separation.size() == 1)
    {
        // Single value for all.
        forAllIter(Map<vector>, field, iter)
        {
            iter() += separation[0];
        }
    }
    else if (separation.size() == field.size())
    {
        forAllIter(Map<vector>, field, iter)
        {
            iter() += separation[iter.key()];
        }
    }
    else
    {
        FatalErrorInFunction
            << "Sizes of field and transformation not equal. field:"
            << field.size() << " transformation:" << separation.size()
            << abort(FatalError);
    }
}


template <>
void Foam::foamSyncTools::separateList
(
    const vectorField& separation,
    EdgeMap<vector>& field
)
{
    if (separation.size() == 1)
    {
        // Single value for all.
        forAllIter(EdgeMap<vector>, field, iter)
        {
            iter() += separation[0];
        }
    }
    else
    {
        FatalErrorInFunction
            << "Multiple separation vectors not supported. field:"
            << field.size() << " transformation:" << separation.size()
            << abort(FatalError);
    }
}*/

// ************************************************************************* //
