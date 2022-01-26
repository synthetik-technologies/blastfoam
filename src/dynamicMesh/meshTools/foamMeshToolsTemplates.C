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
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "globalMeshData.H"
#include "contiguous.H"
#include "transform.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class T>
void Foam::foamSyncTools::separateList
(
    const vectorField& separation,
    UList<T>& field
)
{}


template <class T>
void Foam::foamSyncTools::separateList
(
    const vectorField& separation,
    Map<T>& field
)
{}


template <class T>
void Foam::foamSyncTools::separateList
(
    const vectorField& separation,
    EdgeMap<T>& field
)
{}


// Combine val with existing value at index
template <class T, class CombineOp>
void Foam::foamSyncTools::combine
(
    Map<T>& pointValues,
    const CombineOp& cop,
    const label index,
    const T& val
)
{
    typename Map<T>::iterator iter = pointValues.find(index);

    if (iter != pointValues.end())
    {
        cop(iter(), val);
    }
    else
    {
        pointValues.insert(index, val);
    }
}


// Combine val with existing value at (implicit index) e.
template <class T, class CombineOp>
void Foam::foamSyncTools::combine
(
    EdgeMap<T>& edgeValues,
    const CombineOp& cop,
    const edge& index,
    const T& val
)
{
    typename EdgeMap<T>::iterator iter = edgeValues.find(index);

    if (iter != edgeValues.end())
    {
        cop(iter(), val);
    }
    else
    {
        edgeValues.insert(index, val);
    }
}


template <class T, class CombineOp>
void Foam::foamSyncTools::syncPointMap
(
    const polyMesh& mesh,
    Map<T>& pointValues,        // from mesh point label to value
    const CombineOp& cop,
    const bool applySeparation
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (!hasCouples(patches))
    {
        return;
    }

    // Is there any coupled patch with transformation?
    bool hasTransformation = false;

    if (Pstream::parRun())
    {
        // Send

        forAll (patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].nPoints() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                // Get data per patchPoint in neighbouring point numbers.

                const labelList& meshPts = procPatch.meshPoints();
                const labelList& nbrPts = procPatch.nbrPoints();

                // Extract local values. Create map from nbrPoint to value.
                // Note: how small initial size?
                Map<T> patchInfo(meshPts.size() / 20);

                forAll (meshPts, i)
                {
                    typename Map<T>::const_iterator iter =
                        pointValues.find(meshPts[i]);

                    if (iter != pointValues.end())
                    {
                        if (nbrPts[i] >= 0)
                        {
                            patchInfo.insert(nbrPts[i], iter());
                        }
                    }
                }

                OPstream toNeighb
                (
                    Pstream::commsTypes::blocking,
                    procPatch.neighbProcNo()
                );
                toNeighb << patchInfo;
            }
        }


        // Receive and combine.

        forAll (patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].nPoints() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);
                checkTransform(procPatch, applySeparation);

                IPstream fromNb
                (
                    Pstream::commsTypes::blocking,
                    procPatch.neighbProcNo()
                );
                Map<T> nbrPatchInfo(fromNb);

                if (!procPatch.coupled())
                {
                    hasTransformation = true;
                    transformList(procPatch.forwardT(), nbrPatchInfo);
                }
                else if (applySeparation && procPatch.separated())
                {
                    hasTransformation = true;
                    separateList(-procPatch.separation(), nbrPatchInfo);
                }

                const labelList& meshPts = procPatch.meshPoints();

                // Only update those values which come from neighbour

                forAllConstIter
                (
                    typename Map<T>,
                    nbrPatchInfo,
                    nbrIter
                )
                {
                    combine
                    (
                        pointValues,
                        cop,
                        meshPts[nbrIter.key()],
                        nbrIter()
                    );
                }
            }
        }
    }

    // Do the cyclics.
    forAll (patches, patchI)
    {
        if (isA<cyclicPolyPatch>(patches[patchI]))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchI]);
            checkTransform(cycPatch, applySeparation);

            const edgeList& coupledPoints = cycPatch.coupledPoints();
            const labelList& meshPts = cycPatch.meshPoints();

            // Extract local values. Create map from nbrPoint to value.
            Map<T> half0Values(meshPts.size() / 20);
            Map<T> half1Values(meshPts.size() / 20);

            forAll (coupledPoints, i)
            {
                const edge& e = coupledPoints[i];

                typename Map<T>::const_iterator point0Fnd =
                    pointValues.find(meshPts[e[0]]);

                if (point0Fnd != pointValues.end())
                {
                    half0Values.insert(i, point0Fnd());
                }

                typename Map<T>::const_iterator point1Fnd =
                    pointValues.find(meshPts[e[1]]);

                if (point1Fnd != pointValues.end())
                {
                    half1Values.insert(i, point1Fnd());
                }
            }

            if (!cycPatch.parallel())
            {
                hasTransformation = true;
                transformList(cycPatch.reverseT(), half0Values);
                transformList(cycPatch.forwardT(), half1Values);
            }
            else if (applySeparation && cycPatch.separated())
            {
                hasTransformation = true;

                const vectorField& v = cycPatch.coupledPolyPatch::separation();
                separateList(v, half0Values);
                separateList(-v, half1Values);
            }

            forAll (coupledPoints, i)
            {
                const edge& e = coupledPoints[i];

                typename Map<T>::const_iterator half1Fnd = half1Values.find(i);

                if (half1Fnd != half1Values.end())
                {
                    combine
                    (
                        pointValues,
                        cop,
                        meshPts[e[0]],
                        half1Fnd()
                    );
                }

                typename Map<T>::const_iterator half0Fnd = half0Values.find(i);

                if (half0Fnd != half0Values.end())
                {
                    combine
                    (
                        pointValues,
                        cop,
                        meshPts[e[1]],
                        half0Fnd()
                    );
                }
            }
        }
    }

    // Note: hasTransformation is only used for warning messages so
    // reduction not strictly nessecary.
    // reduce(hasTransformation, orOp<bool>());

    // Synchronize multiple shared points.
    const globalMeshData& pd = mesh.globalData();

    if (pd.nGlobalPoints() > 0)
    {
        if (hasTransformation)
        {
            WarningInFunction
                << "There are decomposed cyclics in this mesh with"
                << " transformations." << endl
                << "This is not supported. The result will be incorrect"
                << endl;
        }
        // meshPoint per local index
        const labelList& sharedPtLabels = pd.sharedPointLabels();
        // global shared index per local index
        const labelList& sharedPtAddr = pd.sharedPointAddr();

        // Values on shared points. Keyed on global shared index.
        Map<T> sharedPointValues(sharedPtAddr.size());


        // Fill my entries in the shared points
        forAll (sharedPtLabels, i)
        {
            label meshPointI = sharedPtLabels[i];

            typename Map<T>::const_iterator fnd =
                pointValues.find(meshPointI);

            if (fnd != pointValues.end())
            {
                combine
                (
                    sharedPointValues,
                    cop,
                    sharedPtAddr[i],    // index
                    fnd()               // value
                );
            }
        }

        // Reduce on master.

        if (Pstream::parRun())
        {
            if (Pstream::master())
            {
                // Receive the edges using shared points from the slave.
                for
                (
                    int slave=Pstream::firstSlave();
                    slave<=Pstream::lastSlave();
                    slave++
                )
                {
                    IPstream fromSlave(Pstream::blocking, slave);
                    Map<T> nbrValues(fromSlave);

                    // Merge neighbouring values with my values
                    forAllConstIter(typename Map<T>, nbrValues, iter)
                    {
                        combine
                        (
                            sharedPointValues,
                            cop,
                            iter.key(), // edge
                            iter()      // value
                        );
                    }
                }

                // Send back
                for
                (
                    int slave=Pstream::firstSlave();
                    slave<=Pstream::lastSlave();
                    slave++
                )
                {
                    OPstream toSlave(Pstream::blocking, slave);
                    toSlave << sharedPointValues;
                }
            }
            else
            {
                // Send to master
                {
                    OPstream toMaster
                    (
                        Pstream::blocking,
                        Pstream::masterNo()
                    );
                    toMaster << sharedPointValues;
                }
                // Receive merged values
                {
                    IPstream fromMaster
                    (
                        Pstream::blocking,
                        Pstream::masterNo()
                    );
                    fromMaster >> sharedPointValues;
                }
            }
        }


        // Merge sharedPointValues (keyed on sharedPointAddr) into
        // pointValues (keyed on mesh points).

        // Map from global shared index to meshpoint
        Map<label> sharedToMeshPoint(2*sharedPtAddr.size());
        forAll (sharedPtAddr, i)
        {
            sharedToMeshPoint.insert(sharedPtAddr[i], sharedPtLabels[i]);
        }

        forAllConstIter(Map<label>, sharedToMeshPoint, iter)
        {
            // Do I have a value for my shared point
            typename Map<T>::const_iterator sharedFnd =
                sharedPointValues.find(iter.key());

            if (sharedFnd != sharedPointValues.end())
            {
                combine
                (
                    pointValues,
                    cop,
                    iter(),     // index
                    sharedFnd() // value
                );
            }
        }
    }
}


template <class T, class CombineOp>
void Foam::foamSyncTools::syncEdgeMap
(
    const polyMesh& mesh,
    EdgeMap<T>& edgeValues,
    const CombineOp& cop,
    const bool applySeparation
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (!hasCouples(patches))
    {
        return;
    }


    // Do synchronisation without constructing globalEdge addressing
    // (since this constructs mesh edge addressing)


    // Swap proc patch info
    // ~~~~~~~~~~~~~~~~~~~~

    if (Pstream::parRun())
    {
        // Send

        forAll (patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].nEdges() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                // Get data per patch edge in neighbouring edge.

                const edgeList& edges = procPatch.edges();
                const labelList& meshPts = procPatch.meshPoints();
                const labelList& nbrPts = procPatch.neighbPoints();

                EdgeMap<T> patchInfo(edges.size() / 20);

                forAll (edges, i)
                {
                    const edge& e = edges[i];
                    const edge meshEdge(meshPts[e[0]], meshPts[e[1]]);

                    typename EdgeMap<T>::const_iterator iter =
                        edgeValues.find(meshEdge);

                    if (iter != edgeValues.end())
                    {
                        const edge nbrEdge(nbrPts[e[0]], nbrPts[e[1]]);

                        if (nbrEdge[0] >= 0 && nbrEdge[1] >= 0)
                        {
                            patchInfo.insert(nbrEdge, iter());
                        }
                    }
                }

                OPstream toNeighb(Pstream::blocking, procPatch.neighbProcNo());
                toNeighb << patchInfo;
            }
        }


        // Receive and combine.

        forAll (patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].nEdges() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);
                checkTransform(procPatch, applySeparation);

                const labelList& meshPts = procPatch.meshPoints();

                IPstream fromNbr(Pstream::blocking, procPatch.neighbProcNo());
                EdgeMap<T> nbrPatchInfo(fromNbr);

                if (!procPatch.parallel())
                {
                    transformList(procPatch.forwardT(), nbrPatchInfo);
                }
                else if (applySeparation && procPatch.separated())
                {
                    separateList(-procPatch.separation(), nbrPatchInfo);
                }

                // Only update those values which come from neighbour

                forAllConstIter
                (
                    typename EdgeMap<T>,
                    nbrPatchInfo,
                    nbrIter
                )
                {
                    const edge& e = nbrIter.key();
                    const edge meshEdge(meshPts[e[0]], meshPts[e[1]]);

                    combine
                    (
                        edgeValues,
                        cop,
                        meshEdge,   // edge
                        nbrIter()   // value
                    );
                }
            }
        }
    }


    // Swap cyclic info
    // ~~~~~~~~~~~~~~~~

    forAll (patches, patchI)
    {
        if (isA<cyclicPolyPatch>(patches[patchI]))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchI]);
            checkTransform(cycPatch, applySeparation);

            const edgeList& coupledEdges = cycPatch.coupledEdges();
            const labelList& meshPts = cycPatch.meshPoints();
            const edgeList& edges = cycPatch.edges();

            // Extract local values. Create map from nbrPoint to value.
            Map<T> half0Values(meshPts.size() / 20);
            Map<T> half1Values(meshPts.size() / 20);

            forAll (coupledEdges, i)
            {
                const edge& twoEdges = coupledEdges[i];

                {
                    const edge& e0 = edges[twoEdges[0]];
                    const edge meshEdge0(meshPts[e0[0]], meshPts[e0[1]]);

                    typename EdgeMap<T>::const_iterator iter =
                        edgeValues.find(meshEdge0);

                    if (iter != edgeValues.end())
                    {
                        half0Values.insert(i, iter());
                    }
                }
                {
                    const edge& e1 = edges[twoEdges[1]];
                    const edge meshEdge1(meshPts[e1[0]], meshPts[e1[1]]);

                    typename EdgeMap<T>::const_iterator iter =
                        edgeValues.find(meshEdge1);

                    if (iter != edgeValues.end())
                    {
                        half1Values.insert(i, iter());
                    }
                }
            }


            // Transform

            if (!cycPatch.parallel())
            {
                transformList(cycPatch.reverseT(), half0Values);
                transformList(cycPatch.forwardT(), half1Values);
            }
            else if (applySeparation && cycPatch.separated())
            {
                const vectorField& v = cycPatch.coupledPolyPatch::separation();
                separateList(v, half0Values);
                separateList(-v, half1Values);
            }


            // Extract and combine information

            forAll (coupledEdges, i)
            {
                const edge& twoEdges = coupledEdges[i];

                typename Map<T>::const_iterator half1Fnd =
                    half1Values.find(i);

                if (half1Fnd != half1Values.end())
                {
                    const edge& e0 = edges[twoEdges[0]];
                    const edge meshEdge0(meshPts[e0[0]], meshPts[e0[1]]);

                    combine
                    (
                        edgeValues,
                        cop,
                        meshEdge0,  // edge
                        half1Fnd()  // value
                    );
                }

                typename Map<T>::const_iterator half0Fnd =
                    half0Values.find(i);
                if (half0Fnd != half0Values.end())
                {
                    const edge& e1 = edges[twoEdges[1]];
                    const edge meshEdge1(meshPts[e1[0]], meshPts[e1[1]]);

                    combine
                    (
                        edgeValues,
                        cop,
                        meshEdge1,  // edge
                        half0Fnd()  // value
                    );
                }
            }
        }
    }

    // Synchronize multiple shared points.
    // Problem is that we don't want to construct shared edges so basically
    // we do here like globalMeshData but then using sparse edge representation
    // (EdgeMap instead of mesh.edges())

    const globalMeshData& pd = mesh.globalData();
    const labelList& sharedPtAddr = pd.sharedPointAddr();
    const labelList& sharedPtLabels = pd.sharedPointLabels();

    // 1. Create map from meshPoint to globalShared index.
    Map<label> meshToShared(2*sharedPtLabels.size());
    forAll (sharedPtLabels, i)
    {
        meshToShared.insert(sharedPtLabels[i], sharedPtAddr[i]);
    }

    // Values on shared points. From two sharedPtAddr indices to a value.
    EdgeMap<T> sharedEdgeValues(meshToShared.size());

    // From shared edge to mesh edge. Used for merging later on.
    EdgeMap<edge> potentialSharedEdge(meshToShared.size());

    // 2. Find any edges using two global shared points. These will always be
    // on the outside of the mesh. (though might not be on coupled patch
    // if is single edge and not on coupled face)
    // Store value (if any) on sharedEdgeValues
    for (label faceI = mesh.nInternalFaces(); faceI < mesh.nFaces(); faceI++)
    {
        const face& f = mesh.faces()[faceI];

        forAll (f, fp)
        {
            label v0 = f[fp];
            label v1 = f[f.fcIndex(fp)];

            Map<label>::const_iterator v0Fnd = meshToShared.find(v0);

            if (v0Fnd != meshToShared.end())
            {
                Map<label>::const_iterator v1Fnd = meshToShared.find(v1);

                if (v1Fnd != meshToShared.end())
                {
                    const edge meshEdge(v0, v1);

                    // edge in shared point labels
                    const edge sharedEdge(v0Fnd(), v1Fnd());

                    // Store mesh edge as a potential shared edge.
                    potentialSharedEdge.insert(sharedEdge, meshEdge);

                    typename EdgeMap<T>::const_iterator edgeFnd =
                        edgeValues.find(meshEdge);

                    if (edgeFnd != edgeValues.end())
                    {
                        // edge exists in edgeValues. See if already in map
                        // (since on same processor, e.g. cyclic)
                        combine
                        (
                            sharedEdgeValues,
                            cop,
                            sharedEdge, // edge
                            edgeFnd()   // value
                        );
                    }
                }
            }
        }
    }


    // Now sharedEdgeValues will contain per potential sharedEdge the value.
    // (potential since an edge having two shared points is not nessecary a
    //  shared edge).
    // Reduce this on the master.

    if (Pstream::parRun())
    {
        if (Pstream::master())
        {
            // Receive the edges using shared points from the slave.
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                IPstream fromSlave(Pstream::blocking, slave);
                EdgeMap<T> nbrValues(fromSlave);

                // Merge neighbouring values with my values
                forAllConstIter(typename EdgeMap<T>, nbrValues, iter)
                {
                    combine
                    (
                        sharedEdgeValues,
                        cop,
                        iter.key(), // edge
                        iter()      // value
                    );
                }
            }

            // Send back
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {

                OPstream toSlave(Pstream::blocking, slave);
                toSlave << sharedEdgeValues;
            }
        }
        else
        {
            // Send to master
            {
                OPstream toMaster(Pstream::blocking, Pstream::masterNo());
                toMaster << sharedEdgeValues;
            }
            // Receive merged values
            {
                IPstream fromMaster(Pstream::blocking, Pstream::masterNo());
                fromMaster >> sharedEdgeValues;
            }
        }
    }


    // Merge sharedEdgeValues (keyed on sharedPointAddr) into edgeValues
    // (keyed on mesh points).

    // Loop over all my shared edges.
    forAllConstIter(typename EdgeMap<edge>, potentialSharedEdge, iter)
    {
        const edge& sharedEdge = iter.key();
        const edge& meshEdge = iter();

        // Do I have a value for the shared edge?
        typename EdgeMap<T>::const_iterator sharedFnd =
            sharedEdgeValues.find(sharedEdge);

        if (sharedFnd != sharedEdgeValues.end())
        {
            combine
            (
                edgeValues,
                cop,
                meshEdge,       // edge
                sharedFnd()     // value
            );
        }
    }
}


template <class T, class CombineOp>
void Foam::foamSyncTools::syncPointList
(
    const polyMesh& mesh,
    UList<T>& pointValues,
    const CombineOp& cop,
    const T& nullValue,
    const bool applySeparation
)
{
    if (pointValues.size() != mesh.nPoints())
    {
        FatalErrorInFunction
            << "Number of values " << pointValues.size()
            << " is not equal to the number of points in the mesh "
            << mesh.nPoints() << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (!hasCouples(patches))
    {
        return;
    }

    // Is there any coupled patch with transformation?
    bool hasTransformation = false;

    if (Pstream::parRun())
    {
        // Send

        forAll (patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].nPoints() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                // Get data per patchPoint in neighbouring point numbers.
                List<T> patchInfo(procPatch.nPoints(), nullValue);

                const labelList& meshPts = procPatch.meshPoints();
                const labelList& nbrPts = procPatch.neighbPoints();

                forAll (nbrPts, pointI)
                {
                    label nbrPointI = nbrPts[pointI];
                    if (nbrPointI >= 0 && nbrPointI < patchInfo.size())
                    {
                        patchInfo[nbrPointI] = pointValues[meshPts[pointI]];
                    }
                }

                OPstream toNbr(Pstream::blocking, procPatch.neighbProcNo());
                toNbr << patchInfo;
            }
        }


        // Receive and combine.

        forAll (patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].nPoints() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);
                checkTransform(procPatch, applySeparation);

                List<T> nbrPatchInfo(procPatch.nPoints());
                {
                    // We do not know the number of points on the other side
                    // so cannot use Pstream::read.
                    IPstream fromNbr
                    (
                        Pstream::blocking,
                        procPatch.neighbProcNo()
                    );
                    fromNbr >> nbrPatchInfo;
                }
                // Null any value which is not on neighbouring processor
                nbrPatchInfo.setSize(procPatch.nPoints(), nullValue);

                if (!procPatch.parallel())
                {
                    hasTransformation = true;
                    transformList(procPatch.forwardT(), nbrPatchInfo);
                }
                else if (applySeparation && procPatch.separated())
                {
                    hasTransformation = true;
                    separateList(-procPatch.separation(), nbrPatchInfo);
                }

                const labelList& meshPts = procPatch.meshPoints();

                forAll (meshPts, pointI)
                {
                    label meshPointI = meshPts[pointI];
                    cop(pointValues[meshPointI], nbrPatchInfo[pointI]);
                }
            }
        }
    }

    // Do the cyclics.
    forAll (patches, patchI)
    {
        if (isA<cyclicPolyPatch>(patches[patchI]))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchI]);

            checkTransform(cycPatch, applySeparation);

            const edgeList& coupledPoints = cycPatch.coupledPoints();
            const labelList& meshPts = cycPatch.meshPoints();

            List<T> half0Values(coupledPoints.size());
            List<T> half1Values(coupledPoints.size());

            forAll (coupledPoints, i)
            {
                const edge& e = coupledPoints[i];

                label point0 = meshPts[e[0]];
                label point1 = meshPts[e[1]];

                half0Values[i] = pointValues[point0];
                half1Values[i] = pointValues[point1];
            }

            if (!cycPatch.parallel())
            {
                hasTransformation = true;
                transformList(cycPatch.reverseT(), half0Values);
                transformList(cycPatch.forwardT(), half1Values);
            }
            else if (applySeparation && cycPatch.separated())
            {
                hasTransformation = true;
                const vectorField& v = cycPatch.coupledPolyPatch::separation();
                separateList(v, half0Values);
                separateList(-v, half1Values);
            }

            forAll (coupledPoints, i)
            {
                const edge& e = coupledPoints[i];

                label point0 = meshPts[e[0]];
                label point1 = meshPts[e[1]];

                cop(pointValues[point0], half1Values[i]);
                cop(pointValues[point1], half0Values[i]);
            }
        }
    }

    //- Note: hasTransformation is only used for warning messages so
    //  reduction not strictly nessecary.
    //reduce(hasTransformation, orOp<bool>());

    // Synchronize multiple shared points.
    const globalMeshData& pd = mesh.globalData();

    if (pd.nGlobalPoints() > 0)
    {
        if (hasTransformation)
        {
            WarningInFunction
                << "There are decomposed cyclics in this mesh with"
                << " transformations." << endl
                << "This is not supported. The result will be incorrect"
                << endl;
        }


        // Values on shared points.
        List<T> sharedPts(pd.nGlobalPoints(), nullValue);

        forAll (pd.sharedPointLabels(), i)
        {
            label meshPointI = pd.sharedPointLabels()[i];
            // Fill my entries in the shared points
            sharedPts[pd.sharedPointAddr()[i]] = pointValues[meshPointI];
        }

        // Combine on master.
        Pstream::listCombineGather(sharedPts, cop);
        Pstream::listCombineScatter(sharedPts);

        // Now we will all have the same information. Merge it back with
        // my local information.
        forAll (pd.sharedPointLabels(), i)
        {
            label meshPointI = pd.sharedPointLabels()[i];
            pointValues[meshPointI] = sharedPts[pd.sharedPointAddr()[i]];
        }
    }
}


template <class T, class CombineOp>
void Foam::foamSyncTools::syncPointList
(
    const polyMesh& mesh,
    const labelList& meshPoints,
    UList<T>& pointValues,
    const CombineOp& cop,
    const T& nullValue,
    const bool applySeparation
)
{
    if (pointValues.size() != meshPoints.size())
    {
        FatalErrorInFunction
            << "Number of values " << pointValues.size()
            << " is not equal to the number of points "
            << meshPoints.size() << abort(FatalError);
    }

    if (!hasCouples(mesh.boundaryMesh()))
    {
        return;
    }

    List<T> meshValues(mesh.nPoints(), nullValue);

    forAll (meshPoints, i)
    {
        meshValues[meshPoints[i]] = pointValues[i];
    }

    foamSyncTools::syncPointList
    (
        mesh,
        meshValues,
        cop,            // combine op
        nullValue,      // null value
        applySeparation // separation
    );

    forAll (meshPoints, i)
    {
        pointValues[i] = meshValues[meshPoints[i]];
    }
}


template <class T, class CombineOp>
void Foam::foamSyncTools::syncEdgeList
(
    const polyMesh& mesh,
    UList<T>& edgeValues,
    const CombineOp& cop,
    const T& nullValue,
    const bool applySeparation
)
{
    if (edgeValues.size() != mesh.nEdges())
    {
        FatalErrorInFunction
            << "Number of values " << edgeValues.size()
            << " is not equal to the number of edges in the mesh "
            << mesh.nEdges() << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (!hasCouples(patches))
    {
        return;
    }

    // Is there any coupled patch with transformation?
    bool hasTransformation = false;

    if (Pstream::parRun())
    {
        // Send

        forAll (patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].nEdges() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                const labelList& meshEdges = procPatch.meshEdges();
                const labelList& neighbEdges = procPatch.neighbEdges();

                // Get region per patch edge in neighbouring edge numbers.
                List<T> patchInfo(procPatch.nEdges(), nullValue);

                forAll (neighbEdges, edgeI)
                {
                    label nbrEdgeI = neighbEdges[edgeI];

                    if (nbrEdgeI >= 0 && nbrEdgeI < patchInfo.size())
                    {
                        patchInfo[nbrEdgeI] = edgeValues[meshEdges[edgeI]];
                    }
                }

                OPstream toNbr(Pstream::blocking, procPatch.neighbProcNo());
                toNbr << patchInfo;
            }
        }

        // Receive and combine.

        forAll (patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].nEdges() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                checkTransform(procPatch, applySeparation);

                const labelList& meshEdges = procPatch.meshEdges();

                // Receive from neighbour. Is per patch edge the region of the
                // neighbouring patch edge.
                List<T> nbrPatchInfo(procPatch.nEdges());

                {
                    IPstream fromNeighb
                    (
                        Pstream::blocking,
                        procPatch.neighbProcNo()
                    );
                    fromNeighb >> nbrPatchInfo;
                }
                // Null any value which is not on neighbouring processor
                nbrPatchInfo.setSize(procPatch.nEdges(), nullValue);

                if (!procPatch.parallel())
                {
                    hasTransformation = true;
                    transformList(procPatch.forwardT(), nbrPatchInfo);
                }
                else if (applySeparation && procPatch.separated())
                {
                    hasTransformation = true;
                    separateList(-procPatch.separation(), nbrPatchInfo);
                }

                forAll (meshEdges, edgeI)
                {
                    label meshEdgeI = meshEdges[edgeI];

                    cop(edgeValues[meshEdgeI], nbrPatchInfo[edgeI]);
                }
            }
        }
    }

    // Do the cyclics.
    forAll (patches, patchI)
    {
        if (isA<cyclicPolyPatch>(patches[patchI]))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchI]);

            checkTransform(cycPatch, applySeparation);

            const edgeList& coupledEdges = cycPatch.coupledEdges();
            const labelList& meshEdges = cycPatch.meshEdges();

            List<T> half0Values(coupledEdges.size());
            List<T> half1Values(coupledEdges.size());

            forAll (coupledEdges, i)
            {
                const edge& e = coupledEdges[i];

                label meshEdge0 = meshEdges[e[0]];
                label meshEdge1 = meshEdges[e[1]];

                half0Values[i] = edgeValues[meshEdge0];
                half1Values[i] = edgeValues[meshEdge1];
            }

            if (!cycPatch.parallel())
            {
                hasTransformation = true;
                transformList(cycPatch.reverseT(), half0Values);
                transformList(cycPatch.forwardT(), half1Values);
            }
            else if (applySeparation && cycPatch.separated())
            {
                hasTransformation = true;

                const vectorField& v = cycPatch.coupledPolyPatch::separation();
                separateList(v, half0Values);
                separateList(-v, half1Values);
            }

            forAll (coupledEdges, i)
            {
                const edge& e = coupledEdges[i];

                label meshEdge0 = meshEdges[e[0]];
                label meshEdge1 = meshEdges[e[1]];

                cop(edgeValues[meshEdge0], half1Values[i]);
                cop(edgeValues[meshEdge1], half0Values[i]);
            }
        }
    }

    //- Note: hasTransformation is only used for warning messages so
    //  reduction not strictly nessecary.
    //reduce(hasTransformation, orOp<bool>());

    // Do the multiple shared edges
    const globalMeshData& pd = mesh.globalData();

    if (pd.nGlobalEdges() > 0)
    {
        if (hasTransformation)
        {
            WarningInFunction
                << "There are decomposed cyclics in this mesh with"
                << " transformations." << endl
                << "This is not supported. The result will be incorrect"
                << endl;
        }

        // Values on shared edges.
        List<T> sharedPts(pd.nGlobalEdges(), nullValue);

        forAll (pd.sharedEdgeLabels(), i)
        {
            label meshEdgeI = pd.sharedEdgeLabels()[i];

            // Fill my entries in the shared edges
            sharedPts[pd.sharedEdgeAddr()[i]] = edgeValues[meshEdgeI];
        }

        // Combine on master.
        Pstream::listCombineGather(sharedPts, cop);
        Pstream::listCombineScatter(sharedPts);

        // Now we will all have the same information. Merge it back with
        // my local information.
        forAll (pd.sharedEdgeLabels(), i)
        {
            label meshEdgeI = pd.sharedEdgeLabels()[i];
            edgeValues[meshEdgeI] = sharedPts[pd.sharedEdgeAddr()[i]];
        }
    }
}


template <class T, class CombineOp>
void Foam::foamSyncTools::syncBoundaryFaceList
(
    const polyMesh& mesh,
    UList<T>& faceValues,
    const CombineOp& cop,
    const bool applySeparation
)
{
    const label nBFaces = mesh.nFaces() - mesh.nInternalFaces();

    if (faceValues.size() != nBFaces)
    {
        FatalErrorInFunction
            << "Number of values " << faceValues.size()
            << " is not equal to the number of boundary faces in the mesh "
            << nBFaces << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (!hasCouples(patches))
    {
        return;
    }


    if (Pstream::parRun())
    {
        // Send

        forAll (patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].size() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                label patchStart = procPatch.start() - mesh.nInternalFaces();

                if (contiguous<T>())
                {
                    OPstream::write
                    (
                        Pstream::blocking,
                        procPatch.neighbProcNo(),
                        reinterpret_cast<const char*>(&faceValues[patchStart]),
                        procPatch.size()*sizeof(T)
                    );
                }
                else
                {
                    OPstream toNbr
                    (
                        Pstream::blocking,
                        procPatch.neighbProcNo()
                    );

                    toNbr <<
                        SubList<T>(faceValues, procPatch.size(), patchStart);
                }
            }
        }


        // Receive and combine.

        forAll (patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].size() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                List<T> nbrPatchInfo(procPatch.size());

                if (contiguous<T>())
                {
                    IPstream::read
                    (
                        Pstream::blocking,
                        procPatch.neighbProcNo(),
                        reinterpret_cast<char*>(nbrPatchInfo.begin()),
                        nbrPatchInfo.byteSize()
                    );
                }
                else
                {
                    IPstream fromNeighb
                    (
                        Pstream::blocking,
                        procPatch.neighbProcNo()
                    );
                    fromNeighb >> nbrPatchInfo;
                }

                if (!procPatch.parallel())
                {
                    transformList(procPatch.forwardT(), nbrPatchInfo);
                }
                else if (applySeparation && procPatch.separated())
                {
                    separateList(-procPatch.separation(), nbrPatchInfo);
                }


                label bFaceI = procPatch.start()-mesh.nInternalFaces();

                forAll (nbrPatchInfo, i)
                {
                    cop(faceValues[bFaceI++], nbrPatchInfo[i]);
                }
            }
        }
    }

    // Do the cyclics.
    forAll (patches, patchI)
    {
        if (isA<cyclicPolyPatch>(patches[patchI]))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchI]);

            label patchStart = cycPatch.start()-mesh.nInternalFaces();

            label half = cycPatch.size()/2;
            label half1Start = patchStart+half;

            List<T> half0Values(SubList<T>(faceValues, half, patchStart));
            List<T> half1Values(SubList<T>(faceValues, half, half1Start));

            if (!cycPatch.parallel())
            {
                transformList(cycPatch.reverseT(), half0Values);
                transformList(cycPatch.forwardT(), half1Values);
            }
            else if (applySeparation && cycPatch.separated())
            {
                const vectorField& v = cycPatch.coupledPolyPatch::separation();
                separateList(v, half0Values);
                separateList(-v, half1Values);
            }

            label i0 = patchStart;
            forAll (half1Values, i)
            {
                cop(faceValues[i0++], half1Values[i]);
            }

            label i1 = half1Start;
            forAll (half0Values, i)
            {
                cop(faceValues[i1++], half0Values[i]);
            }
        }
    }
}


template <class T, class CombineOp>
void Foam::foamSyncTools::syncFaceList
(
    const polyMesh& mesh,
    UList<T>& faceValues,
    const CombineOp& cop,
    const bool applySeparation
)
{
    if (faceValues.size() != mesh.nFaces())
    {
        FatalErrorInFunction
            << "Number of values " << faceValues.size()
            << " is not equal to the number of faces in the mesh "
            << mesh.nFaces() << abort(FatalError);
    }

    SubList<T> bndValues
    (
        faceValues,
        mesh.nFaces()-mesh.nInternalFaces(),
        mesh.nInternalFaces()
    );

    syncBoundaryFaceList
    (
        mesh,
        bndValues,
        cop,
        applySeparation
    );
}


template <class T>
void Foam::foamSyncTools::swapBoundaryFaceList
(
    const polyMesh& mesh,
    UList<T>& faceValues,
    const bool applySeparation
)
{
    syncBoundaryFaceList(mesh, faceValues, eqOp<T>(), applySeparation);
}


template <class T>
void Foam::foamSyncTools::swapFaceList
(
    const polyMesh& mesh,
    UList<T>& faceValues,
    const bool applySeparation
)
{
    syncFaceList(mesh, faceValues, eqOp<T>(), applySeparation);
}


template <unsigned nBits, class CombineOp>
void Foam::foamSyncTools::syncFaceList
(
    const polyMesh& mesh,
    PackedList<nBits>& faceValues,
    const CombineOp& cop
)
{
    if (faceValues.size() != mesh.nFaces())
    {
        FatalErrorInFunction
            << "Number of values " << faceValues.size()
            << " is not equal to the number of faces in the mesh "
            << mesh.nFaces() << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (!hasCouples(patches))
    {
        return;
    }

    // Patch data (proc patches only).
    List<List<unsigned int> > patchValues(patches.size());

    if (Pstream::parRun())
    {
        // Send

        forAll (patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].size() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                patchValues[patchI].setSize(procPatch.size());
                forAll (procPatch, i)
                {
                    patchValues[patchI][i] =
                        faceValues.get(procPatch.start()+i);
                }

                OPstream toNbr(Pstream::blocking, procPatch.neighbProcNo());
                toNbr << patchValues[patchI];
            }
        }


        // Receive and combine.

        forAll (patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].size() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                {
                    IPstream fromNbr
                    (
                        Pstream::blocking,
                        procPatch.neighbProcNo()
                    );
                    fromNbr >> patchValues[patchI];
                }

                // Combine (bitwise)
                forAll (procPatch, i)
                {
                    unsigned int patchVal = patchValues[patchI][i];
                    label meshFaceI = procPatch.start()+i;
                    unsigned int faceVal = faceValues.get(meshFaceI);
                    cop(faceVal, patchVal);
                    faceValues.set(meshFaceI, faceVal);
                }
            }
        }
    }

    // Do the cyclics.
    forAll (patches, patchI)
    {
        if (isA<cyclicPolyPatch>(patches[patchI]))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchI]);

            label half = cycPatch.size()/2;

            for (label i = 0; i < half; i++)
            {
                label meshFace0 = cycPatch.start()+i;
                unsigned int val0 = faceValues.get(meshFace0);
                label meshFace1 = meshFace0 + half;
                unsigned int val1 = faceValues.get(meshFace1);

                unsigned int t = val0;
                cop(t, val1);
                faceValues.set(meshFace0, t);

                cop(val1, val0);
                faceValues.set(meshFace1, val1);
            }
        }
    }
}


template <unsigned nBits>
void Foam::foamSyncTools::swapFaceList
(
    const polyMesh& mesh,
    PackedList<nBits>& faceValues
)
{
    syncFaceList(mesh, faceValues, eqOp<unsigned int>());
}


template <unsigned nBits, class CombineOp>
void Foam::foamSyncTools::syncPointList
(
    const polyMesh& mesh,
    PackedList<nBits>& pointValues,
    const CombineOp& cop,
    const unsigned int nullValue
)
{
    if (pointValues.size() != mesh.nPoints())
    {
        FatalErrorInFunction
            << "Number of values " << pointValues.size()
            << " is not equal to the number of points in the mesh "
            << mesh.nPoints() << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (!hasCouples(patches))
    {
        return;
    }

    // Patch data (proc patches only).
    List<List<unsigned int> > patchValues(patches.size());

    if (Pstream::parRun())
    {
        // Send

        forAll (patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].nPoints() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                patchValues[patchI].setSize(procPatch.nPoints());
                patchValues[patchI] = nullValue;

                const labelList& meshPts = procPatch.meshPoints();
                const labelList& nbrPts = procPatch.neighbPoints();

                forAll (nbrPts, pointI)
                {
                    label nbrPointI = nbrPts[pointI];
                    if (nbrPointI >= 0 && nbrPointI < procPatch.nPoints())
                    {
                        patchValues[patchI][nbrPointI] =
                            pointValues.get(meshPts[pointI]);
                    }
                }

                OPstream toNbr(Pstream::blocking, procPatch.neighbProcNo());
                toNbr << patchValues[patchI];
            }
        }


        // Receive and combine.

        forAll (patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].nPoints() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                {
                    // We do not know the number of points on the other side
                    // so cannot use Pstream::read.
                    IPstream fromNbr
                    (
                        Pstream::blocking,
                        procPatch.neighbProcNo()
                    );
                    fromNbr >> patchValues[patchI];
                }

                // Null any value which is not on neighbouring processor
                patchValues[patchI].setSize(procPatch.nPoints(), nullValue);

                const labelList& meshPts = procPatch.meshPoints();

                forAll (meshPts, pointI)
                {
                    label meshPointI = meshPts[pointI];
                    unsigned int pointVal = pointValues.get(meshPointI);
                    cop(pointVal, patchValues[patchI][pointI]);
                    pointValues.set(meshPointI, pointVal);
                }
            }
        }
    }

    // Do the cyclics.
    forAll (patches, patchI)
    {
        if (isA<cyclicPolyPatch>(patches[patchI]))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchI]);

            const edgeList& coupledPoints = cycPatch.coupledPoints();
            const labelList& meshPts = cycPatch.meshPoints();

            forAll (coupledPoints, i)
            {
                const edge& e = coupledPoints[i];

                label point0 = meshPts[e[0]];
                label point1 = meshPts[e[1]];

                unsigned int val0 = pointValues.get(point0);
                unsigned int t = val0;
                unsigned int val1 = pointValues.get(point1);

                cop(t, val1);
                pointValues.set(point0, t);
                cop(val1, val0);
                pointValues.set(point1, val1);
            }
        }
    }

    // Synchronize multiple shared points.
    const globalMeshData& pd = mesh.globalData();

    if (pd.nGlobalPoints() > 0)
    {
        // Values on shared points. Use unpacked storage for ease!
        List<unsigned int> sharedPts(pd.nGlobalPoints(), nullValue);

        forAll (pd.sharedPointLabels(), i)
        {
            label meshPointI = pd.sharedPointLabels()[i];
            // Fill my entries in the shared points
            sharedPts[pd.sharedPointAddr()[i]] = pointValues.get(meshPointI);
        }

        // Combine on master.
        Pstream::listCombineGather(sharedPts, cop);
        Pstream::listCombineScatter(sharedPts);

        // Now we will all have the same information. Merge it back with
        // my local information.
        forAll (pd.sharedPointLabels(), i)
        {
            label meshPointI = pd.sharedPointLabels()[i];
            pointValues.set(meshPointI, sharedPts[pd.sharedPointAddr()[i]]);
        }
    }
}


template <unsigned nBits, class CombineOp>
void Foam::foamSyncTools::syncEdgeList
(
    const polyMesh& mesh,
    PackedList<nBits>& edgeValues,
    const CombineOp& cop,
    const unsigned int nullValue
)
{
    if (edgeValues.size() != mesh.nEdges())
    {
        FatalErrorInFunction
            << "Number of values " << edgeValues.size()
            << " is not equal to the number of edges in the mesh "
            << mesh.nEdges() << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (!hasCouples(patches))
    {
        return;
    }

    // Patch data (proc patches only).
    List<List<unsigned int> > patchValues(patches.size());

    if (Pstream::parRun())
    {
        // Send

        forAll (patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].nEdges() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                patchValues[patchI].setSize(procPatch.nEdges(), nullValue);

                const labelList& meshEdges = procPatch.meshEdges();
                const labelList& neighbEdges = procPatch.neighbEdges();

                forAll (neighbEdges, edgeI)
                {
                    label nbrEdgeI = neighbEdges[edgeI];
                    if (nbrEdgeI >= 0 && nbrEdgeI < procPatch.nEdges())
                    {
                        patchValues[patchI][nbrEdgeI] =
                            edgeValues.get(meshEdges[edgeI]);
                    }
                }

                OPstream toNbr(Pstream::blocking, procPatch.neighbProcNo());
                toNbr << patchValues[patchI];
            }
        }


        // Receive and combine.

        forAll (patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].nEdges() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                {
                    IPstream fromNeighb
                    (
                        Pstream::blocking,
                        procPatch.neighbProcNo()
                    );
                    fromNeighb >> patchValues[patchI];
                }

                patchValues[patchI].setSize(procPatch.nEdges(), nullValue);

                const labelList& meshEdges = procPatch.meshEdges();

                forAll (meshEdges, edgeI)
                {
                    unsigned int patchVal = patchValues[patchI][edgeI];
                    label meshEdgeI = meshEdges[edgeI];
                    unsigned int edgeVal = edgeValues.get(meshEdgeI);
                    cop(edgeVal, patchVal);
                    edgeValues.set(meshEdgeI, edgeVal);
                }
            }
        }
    }

    // Do the cyclics.
    forAll (patches, patchI)
    {
        if (isA<cyclicPolyPatch>(patches[patchI]))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchI]);

            const edgeList& coupledEdges = cycPatch.coupledEdges();
            const labelList& meshEdges = cycPatch.meshEdges();

            forAll (coupledEdges, i)
            {
                const edge& e = coupledEdges[i];

                label edge0 = meshEdges[e[0]];
                label edge1 = meshEdges[e[1]];

                unsigned int val0 = edgeValues.get(edge0);
                unsigned int t = val0;
                unsigned int val1 = edgeValues.get(edge1);

                cop(t, val1);
                edgeValues.set(edge0, t);
                cop(val1, val0);
                edgeValues.set(edge1, val1);
            }
        }
    }

    // Synchronize multiple shared edges.
    const globalMeshData& pd = mesh.globalData();

    if (pd.nGlobalEdges() > 0)
    {
        // Values on shared edges. Use unpacked storage for ease!
        List<unsigned int> sharedPts(pd.nGlobalEdges(), nullValue);

        forAll (pd.sharedEdgeLabels(), i)
        {
            label meshEdgeI = pd.sharedEdgeLabels()[i];
            // Fill my entries in the shared edges
            sharedPts[pd.sharedEdgeAddr()[i]] = edgeValues.get(meshEdgeI);
        }

        // Combine on master.
        Pstream::listCombineGather(sharedPts, cop);
        Pstream::listCombineScatter(sharedPts);

        // Now we will all have the same information. Merge it back with
        // my local information.
        forAll (pd.sharedEdgeLabels(), i)
        {
            label meshEdgeI = pd.sharedEdgeLabels()[i];
            edgeValues.set(meshEdgeI, sharedPts[pd.sharedEdgeAddr()[i]]);
        }
    }
}


// ************************************************************************* //
