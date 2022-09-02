/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "newLeastSquaresVolPointInterpolation.H"
#include "fvMesh.H"
#include "volFields.H"
#include "pointFields.H"
#include "demandDrivenData.H"
// #include "faMesh.H"
// #include "processorFaPatch.H"
// #include "areaFields.H"
#include "cyclicPolyPatch.H"
#include "cyclicFvPatch.H"
#include "processorPolyPatch.H"
#include "wedgePolyPatch.H"
#include "symmetryPolyPatch.H"
#include "Map.H"
#include "transform.H"
#include "surfaceMesh.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(newLeastSquaresVolPointInterpolation, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void newLeastSquaresVolPointInterpolation::makePointFaces() const
{
    if (debug)
    {
        Info<< "newLeastSquaresVolPointInterpolation::makePointFaces() : "
            << "constructing point boundary faces addressing"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!pointBndFacesPtr_.empty() || !pointProcFacesPtr_.empty())
    {
        FatalErrorIn
        (
            "newLeastSquaresVolPointInterpolation::makePointFaces() const"
        )
            << "point boundary faces adressing already exist"
            << abort(FatalError);
    }

    const pointField& points = mesh().points();
    const labelListList& pointFaces = mesh().pointFaces();
    const labelListList& pointPoints = mesh().pointPoints();

    // Allocate storage for addressing
    pointBndFacesPtr_.set(new labelListList(points.size()));
    labelListList& pointBndFaces = pointBndFacesPtr_();

    // Allocate storage for addressing
    pointCyclicFacesPtr_.set(new labelListList(points.size()));
    labelListList& pointCyclicFaces = pointCyclicFacesPtr_();

    // Allocate storage for addressing
    pointCyclicBndFacesPtr_.set(new labelListList(points.size()));
    labelListList& pointCyclicBndFaces = pointCyclicBndFacesPtr_();


    // Allocate storage for addressing
    pointProcFacesPtr_.set(new labelListList(points.size()));
    labelListList& pointProcFaces = pointProcFacesPtr_();

    forAll(pointBndFaces, pointI)
    {
        const labelList& curPointFaces = pointFaces[pointI];

        labelHashSet bndFaceSet;
        labelHashSet cyclicFaceSet;
        labelHashSet procFaceSet;

        forAll(curPointFaces, faceI)
        {
            label faceID = curPointFaces[faceI];
            label patchID = mesh().boundaryMesh().whichPatch(faceID);

            if (patchID != -1)
            {
                if
                (
                    mesh().boundaryMesh()[patchID].type()
                 == cyclicPolyPatch::typeName
                )
                {
                    cyclicFaceSet.insert(faceID);
                }
                else if
                (
                    mesh().boundaryMesh()[patchID].type()
                 == processorPolyPatch::typeName
                )
                {
                    FatalErrorIn("makePointFaces()")
                        << "processor patches not allowed!"
                        << abort(FatalError);
                    procFaceSet.insert(faceID);
                }
                else if
                (
                    (
                        mesh().boundaryMesh()[patchID].type()
                     != emptyPolyPatch::typeName
                    )
                 && (
                        mesh().boundaryMesh()[patchID].type()
                     != wedgePolyPatch::typeName
                    )
                )
                {
                    bndFaceSet.insert(faceID);
                }
                else if
                (
                    mesh().boundaryMesh()[patchID].type()
                 == wedgePolyPatch::typeName
                )
                {
                    if (pointAxisEdges().found(pointI))
                    {
                        bndFaceSet.insert(faceID);
                    }
                }
            }
        }

        pointBndFaces[pointI] = bndFaceSet.toc();
        pointCyclicFaces[pointI] = cyclicFaceSet.toc();

        {
            labelList allPointProcFaces = procFaceSet.toc();

            // Check for duplicate proc faces
            vectorField allCentres(allPointProcFaces.size());

            forAll(allCentres, faceI)
            {
                label faceID = allPointProcFaces[faceI];
                label patchID = mesh().boundaryMesh().whichPatch(faceID);
                label start = mesh().boundaryMesh()[patchID].start();
                label localFaceID = faceID - start;

                allCentres[faceI] =
                    mesh().C().boundaryField()[patchID][localFaceID];
            }

            boundBox bb(vectorField(points, pointPoints[pointI]), false);
            scalar tol = 0.001*mag(bb.max() - bb.min());

            vectorField centres(allCentres.size(), vector::zero);
            pointProcFaces[pointI] = labelList(allCentres.size(), -1);

            label nCentres = 0;
            forAll(allCentres, faceI)
            {
                bool duplicate = false;
                for (label i=0; i<nCentres; i++)
                {
                    if
                    (
                        mag
                        (
                            centres[i]
                          - allCentres[faceI]
                        )
                      < tol
                    )
                    {
                        duplicate = true;
                        break;
                    }
                }

                if (!duplicate)
                {
                    centres[nCentres] = allCentres[faceI];
                    pointProcFaces[pointI][nCentres] =
                        allPointProcFaces[faceI];
                    nCentres++;
                }
            }

            pointProcFaces[pointI].setSize(nCentres);
        }
    }

    // Find cyclic boundary faces
    forAll(pointCyclicBndFaces, pointI)
    {
        if (pointCyclicFaces[pointI].size())
        {
            if (pointBndFaces[pointI].size())
            {
                label faceID = pointCyclicFaces[pointI][0];
                label patchID = mesh().boundaryMesh().whichPatch(faceID);

                label start = mesh().boundaryMesh()[patchID].start();
                label localFaceID = faceID - start;

                const cyclicPolyPatch& cycPatch =
                    refCast<const cyclicPolyPatch>
                    (
                        mesh().boundaryMesh()[patchID]
                    );

                const edgeList& coupledPoints = cycPatch.coupledPoints();
                const labelList& meshPoints = cycPatch.meshPoints();

                label sizeby2 = cycPatch.size()/2;

//                 label ngbCycPointI = -1;
                forAll(coupledPoints, pI)
                {
                    if (localFaceID < sizeby2)
                    {
                        if ( pointI == meshPoints[coupledPoints[pI][0]] )
                        {
                            pointCyclicBndFaces[pointI] =
                                pointBndFaces
                                [
                                    meshPoints[coupledPoints[pI][1]]
                                ];
//                             ngbCycPointI = meshPoints[coupledPoints[pI][1]];
                            break;
                        }
                    }
                    else
                    {
                        if ( pointI == meshPoints[coupledPoints[pI][1]] )
                        {
                            pointCyclicBndFaces[pointI] =
                                pointBndFaces
                                [
                                    meshPoints[coupledPoints[pI][0]]
                                ];
//                             ngbCycPointI = meshPoints[coupledPoints[pI][0]];
                            break;
                        }
                    }
                }
            }
        }
    }
}

void newLeastSquaresVolPointInterpolation::makeAxisEdges() const
{
    if (debug)
    {
        Info<< "newLeastSquaresVolPointInterpolation::"
            << "makeAxisEdges() : "
            << "constructing axis edges list"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!axisEdgesPtr_.empty())
    {
        FatalErrorIn
        (
            "newLeastSquaresVolPointInterpolation::"
            "makeAxisEdges() const"
        )
            << "axis edges list already exist"
            << abort(FatalError);
    }

    labelHashSet axisEdgeSet;

    forAll(mesh().boundaryMesh(), patchI)
    {
        if
        (
            mesh().boundaryMesh()[patchI].type()
         == wedgePolyPatch::typeName
        )
        {
            const wedgePolyPatch& wedge =
                refCast<const wedgePolyPatch>(mesh().boundaryMesh()[patchI]);

            const labelList& meshEdges = wedge.meshEdges();

            const labelListList& edgeFaces = mesh().edgeFaces();

            forAll(meshEdges, edgeI)
            {
                if (!wedge.isInternalEdge(edgeI))
                {
                    label curMshEdge = meshEdges[edgeI];

                    const labelList& curEdgeFaces = edgeFaces[curMshEdge];

                    if (curEdgeFaces.size() == 2)
                    {
                        label patch0 =
                            mesh().boundaryMesh().whichPatch
                            (
                                curEdgeFaces[0]
                            );

                        label patch1 =
                            mesh().boundaryMesh().whichPatch
                            (
                                curEdgeFaces[1]
                            );

                        if
                        (
                            mesh().boundaryMesh()[patch0].type()
                         == mesh().boundaryMesh()[patch1].type()
                        )
                        {
                            axisEdgeSet.insert(curMshEdge);
                        }
                    }
                }
            }

            break;
        }
    }

    axisEdgesPtr_.set(new labelList(axisEdgeSet.toc()));
}

void newLeastSquaresVolPointInterpolation::makePointAxisEdges() const
{
    if (debug)
    {
        Info<< "newLeastSquaresVolPointInterpolation::"
            << "makePointAxisEdges() : "
            << "constructing point axis edges addressing"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!pointAxisEdgesPtr_.empty())
    {
        FatalErrorIn
        (
            "newLeastSquaresVolPointInterpolation::"
            "makePointAxisEdges() const"
        )
            << "point axis edges addressing already exists"
            << abort(FatalError);
    }

    pointAxisEdgesPtr_.set(new Map<labelList>());

    Map<labelList>& pointAxisEdges =
        pointAxisEdgesPtr_();

    const edgeList& edges = mesh().edges();

    List<labelHashSet> pae(mesh().points().size());

    forAll(axisEdges(), edgeI)
    {
        label curEdge = axisEdges()[edgeI];

        pae[edges[curEdge].start()].insert(curEdge);
        pae[edges[curEdge].end()].insert(curEdge);
    }

    forAll(pae, pointI)
    {
        labelList curEdges = pae[pointI].toc();

        if (curEdges.size())
        {
            pointAxisEdges.insert(pointI, curEdges);
        }
    }

    if (debug)
    {
        Info<< "point-axis-edges: " << pointAxisEdges << endl;
    }
}

//void newLeastSquaresVolPointInterpolation::makePointNgbProcBndFaceCentres() const
// {
//     if (debug)
//     {
//         Info<< "newLeastSquaresVolPointInterpolation::"
//             << "makePointNgbProcFaceCentres() : "
//             << "constructing point ngb processor face centres"
//             << endl;
//     }

//     // It is an error to attempt to recalculate
//     // if the pointer is already set
//     if (pointNgbProcBndFaceCentresPtr_)
//     {
//         FatalErrorIn
//         (
//             "newLeastSquaresVolPointInterpolation::"
//             "makePointNgbProcFaceCentres() const"
//         )
//             << "point ngb processor face centres already exist"
//             << abort(FatalError);
//     }

//     pointNgbProcBndFaceCentresPtr_ = new Map<vectorField>();

//     Map<vectorField>& pointNgbProcBndFaceCentres =
//         *pointNgbProcBndFaceCentresPtr_;

//     pointNgbProcBndFaceFieldData(mesh().C(), pointNgbProcBndFaceCentres);
// }

void newLeastSquaresVolPointInterpolation::
makeGlobalPointNgbProcBndFaceCentres() const
{
    if (debug)
    {
        Info<< "newLeastSquaresVolPointInterpolation::"
            << "makeGlobalPointNgbProcBndFaceCentres() : "
            << "constructing global point ngb processor bnd face centres"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!globalPointNgbProcBndFaceCentresPtr_.empty())
    {
        FatalErrorIn
        (
            word("newLeastSquaresVolPointInterpolation::")
          + word("makeGlobalPointNgbProcBndFaceCentres() const")
        )
            << "global point ngb processor bnd face centres already exist"
                << abort(FatalError);
    }

    globalPointNgbProcBndFaceCentresPtr_.set(new Map<vectorField>());

    Map<vectorField>& globalPointNgbProcBndFaceCentres =
        globalPointNgbProcBndFaceCentresPtr_();

    globalPointNgbProcBndFaceFieldData
    (
        mesh().C(),
        globalPointNgbProcBndFaceCentres
    );
}

void newLeastSquaresVolPointInterpolation::
makeGlobalPointNgbProcCellCentres() const
{
    if (debug)
    {
        Info<< "newLeastSquaresVolPointInterpolation::"
            << "makeGlobalPointNgbProcCellCentres() : "
            << "constructing global point ngb processor cell centres"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!globalPointNgbProcCellCentresPtr_.empty())
    {
        FatalErrorIn
        (
            word("newLeastSquaresVolPointInterpolation::")
          + word("makeGlobalPointNgbProcCellCentres() const")
        )
            << "global point ngb processor cell centres already exist"
                << abort(FatalError);
    }

    globalPointNgbProcCellCentresPtr_.set(new Map<vectorField>());

    Map<vectorField>& globalPointNgbProcCellCentres =
        globalPointNgbProcCellCentresPtr_();

    globalPointNgbProcCellFieldData(mesh().C(), globalPointNgbProcCellCentres);
}


void newLeastSquaresVolPointInterpolation::makeProcBndFaces() const
{
    if (debug)
    {
        Info<< "newLeastSquaresVolPointInterpolation::makeProcBndFaces() : "
            << "constructing list of boundary faces needed by neighbour "
            << "processors"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!procBndFacesPtr_.empty())
    {
        FatalErrorIn
        (
            "newLeastSquaresVolPointInterpolation::makeProcBndFaces() const"
        )
            << "List of boundary faces needed by neighbour processors "
                << "already exist"
                << abort(FatalError);
    }

    procBndFacesPtr_.set(new labelListList(Pstream::nProcs()));
    labelListList& procBndFaces = procBndFacesPtr_();
    forAll(procBndFaces, procI)
    {
        procBndFaces[procI] = labelList(0);
    }

    pointProcBndFacesPtr_.set
    (
        new List<List<labelPair> >
        (
            mesh().points().size(),
            List<labelPair>(0)
        )
    );
    List<List<labelPair> >& pointProcBndFaces = pointProcBndFacesPtr_();

//     pointProcBndFacesPtr_ = new Map<List<labelPair> >();
//     Map<List<labelPair> >& pointProcBndFaces = *pointProcBndFacesPtr_;

    const labelListList& ptBndFaces = pointBndFaces();

//     const labelListList& pointCells = mesh().pointCells();
//     const vectorField& C = mesh().cellCentres();

    forAll(mesh().boundaryMesh(), patchI)
    {
        if
        (
            mesh().boundaryMesh()[patchI].type()
         == processorPolyPatch::typeName
        )
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>
                (
                    mesh().boundaryMesh()[patchI]
                );

            const labelList& bndPoints = procPatch.boundaryPoints();
            const labelList& meshPoints = procPatch.meshPoints();

            labelHashSet faceSet;
            labelHashSet pointSet;
            SLList<labelPair> pointFaceSet;

            const labelList& glPoints =
                mesh().globalData().sharedPointLabels();

            forAll(bndPoints, pointI)
            {
                label curPoint = bndPoints[pointI];
                label curMeshPoint = meshPoints[curPoint];

                bool glPoint(findIndex(glPoints, curMeshPoint) != -1);

                if (!glPoint)
                {
                    const labelList& curPointBndFaces =
                        ptBndFaces[curMeshPoint];

                    forAll(curPointBndFaces, faceI)
                    {
                        if (!faceSet.found(curPointBndFaces[faceI]))
                        {
                            faceSet.insert(curPointBndFaces[faceI]);
                        }
                        if (!pointSet.found(curPoint))
                        {
                            pointSet.insert(curPoint);
                        }
                        pointFaceSet.insert
                        (
                            labelPair(curPoint, curPointBndFaces[faceI])
                        );
                    }

//                     if (curPointBndFaces.size())
//                     {
//                         Pout<< procPatch.neighbProcNo()
//                             << ", "
//                             << curMeshPoint
//                             << ", "
//                             << pointCells[curMeshPoint]
//                             << ", "
//                             << curPointBndFaces
//                             << ", "
//                             << vectorField(C, pointCells[curMeshPoint]);

//                         forAll(curPointBndFaces, fI)
//                         {
//                             Pout<< ", " <<
//                                 mesh().boundaryMesh().whichPatch
//                                 (
//                                     curPointBndFaces[fI]
//                                 );
//                         }
//                         Pout<< endl;
//                     }
                }
            }

            procBndFaces[procPatch.neighbProcNo()] = faceSet.toc();

            // Point face addressing
            labelList patchPoints = pointSet.toc();

            // List(SLList) in the Foundation version doesn't work, maybe only
            // for Type = labelPair, so we will do it manually
            List<labelPair> patchPointsFaces(pointFaceSet.size());
            {
                label i = 0;
                for
                (
                    typename SLList<labelPair>::const_iterator iter =
                        pointFaceSet.begin();
                    iter != pointFaceSet.end();
                    ++iter
                    )
                {
                    patchPointsFaces[i++] = iter();
                }
            }

            labelListList patchPointFaces(patchPoints.size());

            forAll(patchPoints, pointI)
            {
                labelHashSet faceSet;

                forAll(patchPointsFaces, pI)
                {
                    if
                    (
                        patchPointsFaces[pI].first()
                     == patchPoints[pointI]
                    )
                    {
                        label pointFace =
                            findIndex
                            (
                                procBndFaces[procPatch.neighbProcNo()],
                                patchPointsFaces[pI].second()
                            );
                        faceSet.insert(pointFace);
                    }
                }

                patchPointFaces[pointI] = faceSet.toc();
            }

            // Parallel data exchange
            {
                OPstream toNeighbProc
                (
                    Pstream::commsTypes::blocking,
                    procPatch.neighbProcNo()
                    // size of field
                );

                toNeighbProc << patchPoints << patchPointFaces;
            }

            labelList ngbPatchPoints;
            labelListList ngbPatchPointFaces;

            {
                IPstream fromNeighbProc
                (
                    Pstream::commsTypes::blocking,
                    procPatch.neighbProcNo()
                    // size of field
                );

                fromNeighbProc >> ngbPatchPoints >> ngbPatchPointFaces;
            }

            forAll(ngbPatchPoints, pointI)
            {
                label curNgbPoint = ngbPatchPoints[pointI];

                label curLocalPoint =
                    findIndex(procPatch.neighbPoints(), curNgbPoint);
//                     procPatch.neighbPoints()[curNgbPoint];

                List<labelPair> addressing
                (
                    ngbPatchPointFaces[pointI].size(),
                    labelPair(-1, -1)
                );
                forAll(addressing, faceI)
                {
                    addressing[faceI] =
                        labelPair
                        (
                            procPatch.neighbProcNo(),
                            ngbPatchPointFaces[pointI][faceI]
                        );
                }

                pointProcBndFaces[meshPoints[curLocalPoint]] = addressing;
//                 pointProcBndFaces.insert
//                 (
//                     meshPoints[curLocalPoint],
//                     addressing
//                 );
            }
        }
    }
}


void newLeastSquaresVolPointInterpolation::makeProcBndFaceCentres() const
{
    if (debug)
    {
        Info<< "newLeastSquaresVolPointInterpolation::makeProcBndFaceCentres() : "
            << "constructing centres of boundary faces from ngb processors"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!procBndFaceCentresPtr_.empty())
    {
        FatalErrorIn
        (
            "newLeastSquaresVolPointInterpolation::makeProcBndFaceCentres() const"
        )
            << "Centres of faces from ngb processors already exist"
                << abort(FatalError);
    }

    procBndFaceCentresPtr_.set
    (
        new FieldField<Field, vector>
        (
            Pstream::nProcs()
        )
    );
    FieldField<Field, vector>& procBndFaceCentres = procBndFaceCentresPtr_();

    const vectorField& Cf = mesh().faceCentres();

    if (Pstream::parRun())
    {
        // Send centres
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                vectorField procCentres(Cf, procBndFaces()[procI]);

                {
                    OPstream toNeighbProc
                    (
                        Pstream::commsTypes::blocking,
                        procI
                        // size of field
                    );

                    toNeighbProc << procCentres;
                }
            }
        }

        // Receive centres
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                vectorField procCentres;

                {
                    IPstream fromNeighbProc
                    (
                        Pstream::commsTypes::blocking,
                        procI
                        // size of field
                    );

                    fromNeighbProc >> procCentres;
                }

                procBndFaceCentres.set(procI, new vectorField(procCentres));
            }
            else
            {
                procBndFaceCentres.set(procI, new vectorField(0));
            }
        }
    }
    else
    {
        forAll (procBndFaceCentres, procI)
        {
            procBndFaceCentres.set(procI, new vectorField(0));
        }
    }
}


void newLeastSquaresVolPointInterpolation::makeProcCells() const
{
    if (debug)
    {
        Info<< "newLeastSquaresVolPointInterpolation::makeProcCells() : "
            << "constructing list of cells needed by neighbour processors"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!procCellsPtr_.empty())
    {
        FatalErrorIn
        (
            "newLeastSquaresVolPointInterpolation::makeProcCells() const"
        )
            << "List of cells needed by neighbour processors already exist"
                << abort(FatalError);
    }

    procCellsPtr_.set(new labelListList(Pstream::nProcs()));
    labelListList& procCells = procCellsPtr_();
    forAll(procCells, procI)
    {
        procCells[procI] = labelList(0);
    }

    pointProcCellsPtr_.set(new Map<List<labelPair> >());
    Map<List<labelPair> >& pointProcCells = pointProcCellsPtr_();

    const labelListList& pointCells = mesh().pointCells();

    forAll(mesh().boundaryMesh(), patchI)
    {
        if
        (
            mesh().boundaryMesh()[patchI].type()
         == processorPolyPatch::typeName
        )
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>
                (
                    mesh().boundaryMesh()[patchI]
                );

            const labelList& meshPoints = procPatch.meshPoints();

            const unallocLabelList& patchCells = procPatch.faceCells();

            labelHashSet cellSet;
            labelHashSet pointSet;
            SLList<labelPair> pointCellSet;

            const labelList& glPoints =
                mesh().globalData().sharedPointLabels();

            forAll(meshPoints, pointI)
            {
                label curPoint = meshPoints[pointI];

                bool glPoint(findIndex(glPoints, curPoint) != -1);

                labelHashSet localCellSet;

                if (!glPoint)
                {
                    const labelList& curCells = pointCells[curPoint];

                    forAll(curCells, cellI)
                    {
                        if (findIndex(patchCells, curCells[cellI]) == -1)
                        {
                            if (!cellSet.found(curCells[cellI]))
                            {
                                cellSet.insert(curCells[cellI]);
                            }

                            if (!localCellSet.found(curCells[cellI]))
                            {
                                localCellSet.insert(curCells[cellI]);
                            }

                            if (!pointSet.found(pointI))
                            {
                                pointSet.insert(pointI);
                            }

                            pointCellSet.insert
                            (
                                labelPair(pointI, curCells[cellI])
                            );
                        }
                    }
                }
            }

            procCells[procPatch.neighbProcNo()] = cellSet.toc();

            // Point cell addressing
            labelList patchPoints = pointSet.toc();
            // List(SLList) in the Foundation version doesn't work, maybe only
            // for Type = labelPair, so we will do it manually
            List<labelPair> patchPointsCells(pointCellSet.size());
            {
                label i = 0;
                for
                (
                    typename SLList<labelPair>::const_iterator iter =
                        pointCellSet.begin();
                    iter != pointCellSet.end();
                    ++iter
                    )
                {
                    patchPointsCells[i++] = iter();
                }
            }

            labelListList patchPointCells(patchPoints.size());

            label nCells = 0;

            forAll(patchPointCells, pointI)
            {
                labelHashSet cellSet;

                forAll(patchPointsCells, pI)
                {
                    if
                    (
                        patchPointsCells[pI].first()
                     == patchPoints[pointI]
                    )
                    {
                        label pointCell =
                            findIndex
                            (
                                procCells[procPatch.neighbProcNo()],
                                patchPointsCells[pI].second()
                            );
                        cellSet.insert(pointCell);
                    }
                }

                patchPointCells[pointI] = cellSet.toc();

                nCells += patchPointCells[pointI].size();
            }

            // Parallel data exchange
            {
                OPstream toNeighbProc
                (
                    Pstream::commsTypes::blocking,
                    procPatch.neighbProcNo()
                    // size of field
//                   + 2*patchPoints.size()*sizeof(label)
//                   + nCells*sizeof(label)
                );

                toNeighbProc << patchPoints << patchPointCells;
            }
        }
    }


    forAll(mesh().boundaryMesh(), patchI)
    {
        if
        (
            mesh().boundaryMesh()[patchI].type()
         == processorPolyPatch::typeName
        )
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>
                (
                    mesh().boundaryMesh()[patchI]
                );

            const labelList& meshPoints = procPatch.meshPoints();

            labelList ngbPatchPoints;
            labelListList ngbPatchPointCells;

            {
                IPstream fromNeighbProc
                (
                    Pstream::commsTypes::blocking,
                    procPatch.neighbProcNo()
                    // size of field
                );

                fromNeighbProc >> ngbPatchPoints >> ngbPatchPointCells;
            }

            forAll(ngbPatchPoints, pointI)
            {
                label curNgbPoint = ngbPatchPoints[pointI];

                label curLocalPoint =
                    findIndex(procPatch.neighbPoints(), curNgbPoint);
//                     procPatch.neighbPoints()[curNgbPoint];

                List<labelPair> addressing
                (
                    ngbPatchPointCells[pointI].size(),
                    labelPair(-1, -1)
                );

                forAll(addressing, cellI)
                {
                    addressing[cellI] =
                        labelPair
                        (
                            procPatch.neighbProcNo(),
                            ngbPatchPointCells[pointI][cellI]
                        );
                }

                pointProcCells.insert
                (
                    meshPoints[curLocalPoint],
                    addressing
                );
            }
        }
    }
}


void newLeastSquaresVolPointInterpolation::makeProcCellCentres() const
{
    if (debug)
    {
        Info<< "newLeastSquaresVolPointInterpolation::makeProcCellCentres() : "
            << "constructing centres of cells from ngb processors"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!procCellCentresPtr_.empty())
    {
        FatalErrorIn
        (
            "newLeastSquaresVolPointInterpolation::makeProcCellCentres() const"
        )
            << "Centres of cells from ngb processors already exist"
                << abort(FatalError);
    }

    procCellCentresPtr_.set
    (
        new FieldField<Field, vector>
        (
            Pstream::nProcs()
        )
    );
    FieldField<Field, vector>& procCellCentres = procCellCentresPtr_();

//     const vectorField& CI = mesh().C().internalField();
    const vectorField& CI = mesh().cellCentres();

    if (Pstream::parRun())
    {
        // Send centres
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                vectorField procCentres(CI, procCells()[procI]);

                {
                    OPstream toNeighbProc
                    (
                        Pstream::commsTypes::blocking,
                        procI
                        // size of field
                    );

                    toNeighbProc << procCentres;
                }
            }
        }

        // Receive centres
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                vectorField procCentres;

                {
                    IPstream fromNeighbProc
                    (
                        Pstream::commsTypes::blocking,
                        procI
                        // size of field
                    );

                    fromNeighbProc >> procCentres;
                }

                procCellCentres.set(procI, new vectorField(procCentres));
            }
            else
            {
                procCellCentres.set(procI, new vectorField(0));
            }
        }
    }
    else
    {
        forAll (procCellCentres, procI)
        {
            procCellCentres.set(procI, new vectorField(0));
        }
    }
}


void newLeastSquaresVolPointInterpolation::makeWeights() const
{
    if (debug)
    {
        Info<< "newLeastSquaresVolPointInterpolation::makeWeights() : "
            << "constructing weights"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!weightsPtr_.empty())
    {
        FatalErrorIn
        (
            "newLeastSquaresVolPointInterpolation::makeWeights() const"
        )
            << "weights already exist"
            << abort(FatalError);
    }

    weightsPtr_.set
    (
        new FieldField<Field, scalar>(mesh().points().size())
    );
    FieldField<Field, scalar>& weights = weightsPtr_();

    const vectorField& p = mesh().points();
    const vectorField& C = mesh().cellCentres();
    const vectorField& Cf = mesh().faceCentres();

    const labelListList& ptCells = mesh().pointCells();
    const labelListList& ptBndFaces = pointBndFaces();
    const labelListList& ptCyclicFaces = pointCyclicFaces();
    const labelListList& ptCyclicBndFaces = pointCyclicBndFaces();
    const labelListList& ptProcFaces = pointProcFaces();

//     const Map<vectorField>& ptNgbProcBndFaceCentres =
//         pointNgbProcBndFaceCentres();

    const Map<vectorField>& gPtNgbProcBndFaceCentres =
        globalPointNgbProcBndFaceCentres();

    const Map<vectorField>& gPtNgbProcCellCentres =
        globalPointNgbProcCellCentres();

    const Map<List<labelPair> >& ptProcCells = pointProcCells();

    const List<List<labelPair> >& ptProcBndFaces = pointProcBndFaces();
//     const Map<List<labelPair> >& ptProcBndFaces = pointProcBndFaces();

    const FieldField<Field, vector>& procCentres = procCellCentres();

    const FieldField<Field, vector>& procBndFaceCent = procBndFaceCentres();

    forAll(weights, pointI)
    {
        const labelList& interpCells = ptCells[pointI];
        const labelList& interpBndFaces = ptBndFaces[pointI];
        const labelList& interpCyclicFaces = ptCyclicFaces[pointI];
        const labelList& interpCyclicBndFaces = ptCyclicBndFaces[pointI];
        const labelList& interpProcFaces = ptProcFaces[pointI];

//         vectorField interpNgbProcBndFaceCentres(0);

//         // Boundary faces from neighbour processors
//         if (ptNgbProcBndFaceCentres.found(pointI))
//         {
//              interpNgbProcBndFaceCentres =
//                 ptNgbProcBndFaceCentres[pointI];
//         }

        vectorField glInterpNgbProcBndFaceCentres(0);

        // Boundary faces from neighbour processors
        if (gPtNgbProcBndFaceCentres.found(pointI))
        {
             glInterpNgbProcBndFaceCentres =
                gPtNgbProcBndFaceCentres[pointI];
        }

        vectorField glInterpNgbProcCellCentres(0);

        // Boundary faces from neighbour processors
        if (gPtNgbProcCellCentres.found(pointI))
        {
             glInterpNgbProcCellCentres =
                gPtNgbProcCellCentres[pointI];
        }

        vectorField interpNgbProcCellCentres(0);

        if (ptProcCells.found(pointI))
        {
            const List<labelPair>& pc = ptProcCells[pointI];

            interpNgbProcCellCentres.setSize(pc.size());

            forAll(pc, cI)
            {
                interpNgbProcCellCentres[cI] =
                    procCentres
                    [
                        pc[cI].first()
                    ]
                    [
                        pc[cI].second()
                    ];
            }
        }

        vectorField interpNgbProcBndFaceCentres(0);

        if (ptProcBndFaces[pointI].size())
//         if (ptProcBndFaces.found(pointI))
        {
            const List<labelPair>& pf = ptProcBndFaces[pointI];

            interpNgbProcBndFaceCentres.setSize(pf.size());

            forAll(pf, fI)
            {
                interpNgbProcBndFaceCentres[fI] =
                    procBndFaceCent
                    [
                        pf[fI].first()
                    ]
                    [
                        pf[fI].second()
                    ];
            }
        }

        vectorField allPoints
        (
            interpCells.size()
          + interpBndFaces.size()
          + interpCyclicFaces.size()
          + interpCyclicBndFaces.size()
          + interpProcFaces.size()
//        + interpNgbProcBndFaceCentres.size()
          + glInterpNgbProcBndFaceCentres.size()
          + glInterpNgbProcCellCentres.size()
          + interpNgbProcCellCentres.size()
          + interpNgbProcBndFaceCentres.size(),
            vector::zero
        );

        label pointID = 0;

        // Cells
        for (label i=0; i<interpCells.size(); i++)
        {
            allPoints[pointID++] = C[interpCells[i]];
        }

        // Boundary faces
        for (label i=0; i<interpBndFaces.size(); i++)
        {
            label faceID = interpBndFaces[i];

            allPoints[pointID++] = Cf[faceID];
        }

        // Cyclic boundary faces
        for (label i=0; i<interpCyclicFaces.size(); i++)
        {
            label faceID = interpCyclicFaces[i];
            label patchID = mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            const unallocLabelList& faceCells =
                mesh().boundary()[patchID].faceCells();

            const cyclicFvPatch& cycPatch =
                refCast<const cyclicFvPatch>(mesh().boundary()[patchID]);

            label sizeby2 = faceCells.size()/2;

            if (cycPatch.parallel())
            {
                if (localFaceID < sizeby2)
                {
                    vector deltaNgb = //dni
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID + sizeby2
                        ]
                      - C[faceCells[localFaceID + sizeby2]];

                    vector deltaOwn = //ddi
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID
                        ]
                      - C[faceCells[localFaceID]];

                    vector delta = deltaOwn - deltaNgb;

                    allPoints[pointID++] = C[faceCells[localFaceID]] + delta;
                }
                else
                {
                    vector deltaNgb = //ddi
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID - sizeby2
                        ]
                      - C[faceCells[localFaceID - sizeby2]];

                    vector deltaOwn = //dni
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID
                        ]
                      - C[faceCells[localFaceID]];

                    vector delta = deltaNgb - deltaOwn;
                    delta *= -1;

                    allPoints[pointID++] = C[faceCells[localFaceID]] + delta;
                }
            }
            else
            {
                if (localFaceID < sizeby2)
                {
                    vector deltaNgb = // dni
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID + sizeby2
                        ]
                      - C[faceCells[localFaceID + sizeby2]];

                    vector deltaOwn = //ddi
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID
                        ]
                      - C[faceCells[localFaceID]];

                    vector delta =
                        deltaOwn
                      - transform(cycPatch.forwardT()[0], deltaNgb);

                    allPoints[pointID++] = C[faceCells[localFaceID]] + delta;
                }
                else
                {
                    vector deltaNgb = //ddi
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID - sizeby2
                        ]
                      - C[faceCells[localFaceID - sizeby2]];

                    vector deltaOwn = //dni
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID
                        ]
                      - C[faceCells[localFaceID]];

                    vector delta =
                        deltaNgb
                      - transform(cycPatch.forwardT()[0], deltaOwn);

                    delta = -transform(cycPatch.reverseT()[0], delta);

                    allPoints[pointID++] = C[faceCells[localFaceID]] + delta;
                }
            }
        }

        for (label i=0; i<interpCyclicBndFaces.size(); i++)
        {
            label cycFaceID = interpCyclicFaces[0];
            label cycPatchID = mesh().boundaryMesh().whichPatch(cycFaceID);
            label cycStart = mesh().boundaryMesh()[cycPatchID].start();
            label cycLocalFaceID = cycFaceID - cycStart;

            const cyclicFvPatch& cycPatch =
                refCast<const cyclicFvPatch>(mesh().boundary()[cycPatchID]);

            const cyclicPolyPatch& cycPolyPatch =
                refCast<const cyclicPolyPatch>
                (
                    mesh().boundaryMesh()[cycPatchID]
                );

            label faceID = interpCyclicBndFaces[i];

            const edgeList& coupledPoints = cycPolyPatch.coupledPoints();
            const labelList& meshPoints = cycPolyPatch.meshPoints();

            label sizeby2 = cycPolyPatch.size()/2;

            label ngbCycPointI = -1;
            forAll(coupledPoints, pI)
            {
                if (cycLocalFaceID < sizeby2)
                {
                    if ( pointI == meshPoints[coupledPoints[pI][0]] )
                    {
                        ngbCycPointI = meshPoints[coupledPoints[pI][1]];
                        break;
                    }
                }
                else
                {
                    if ( pointI == meshPoints[coupledPoints[pI][1]] )
                    {
                        ngbCycPointI = meshPoints[coupledPoints[pI][0]];
                        break;
                    }
                }
            }

            if (cycPatch.parallel())
            {
                vector deltaNgb = //dni
                    Cf[faceID] - p[ngbCycPointI];

                allPoints[pointID++] = p[pointI] + deltaNgb;
            }
            else
            {
                vector deltaNgb = //dni
                    Cf[faceID] - p[ngbCycPointI];


                if (cycLocalFaceID < sizeby2)
                {
                    allPoints[pointID++] =
                        p[pointI]
                      + transform(cycPatch.forwardT()[0], deltaNgb);
                }
                else
                {
                    allPoints[pointID++] =
                        p[pointI]
                      + transform(cycPatch.reverseT()[0], deltaNgb);
                }
            }
        }

        // Processor boundary faces
        for (label i=0; i<interpProcFaces.size(); i++)
        {
            label faceID = interpProcFaces[i];
            label patchID = mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            allPoints[pointID++] =
                mesh().C().boundaryField()[patchID][localFaceID];
        }

//         // Boundary faces from neighbour processors
//         for (label i=0; i<interpNgbProcBndFaceCentres.size(); i++)
//         {
//             allPoints[pointID++] =
//                 interpNgbProcBndFaceCentres[i];
//         }

        // Global point bnd faces from neighbour processors
        for (label i=0; i<glInterpNgbProcBndFaceCentres.size(); i++)
        {
            allPoints[pointID++] =
                glInterpNgbProcBndFaceCentres[i];
        }

        // Global point bnd faces from neighbour processors
        for (label i=0; i<glInterpNgbProcCellCentres.size(); i++)
        {
            allPoints[pointID++] =
                glInterpNgbProcCellCentres[i];
        }

        // Cells from neighbour processors
        for (label i=0; i<interpNgbProcCellCentres.size(); i++)
        {
            allPoints[pointID++] =
                interpNgbProcCellCentres[i];
        }

        // Boundary faces from neighbour processors
        for (label i=0; i<interpNgbProcBndFaceCentres.size(); i++)
        {
            allPoints[pointID++] =
                interpNgbProcBndFaceCentres[i];
        }

        vectorField allMirrorPoints(0);
        if (mag(mirrorPlaneTransformation()[pointI].first())>SMALL)
//         if (mirrorPlaneTransformation().found(pointI))
        {
            const vector& n = mirrorPlaneTransformation()[pointI].first();

            allMirrorPoints.setSize(allPoints.size());

            forAll(allPoints, pI)
            {
                allMirrorPoints[pI] =
                    p[pointI] + transform(I-2*n*n, (allPoints[pI]-p[pointI]));
            }
        }

        // Weights
        scalarField W(allPoints.size() + allMirrorPoints.size(), 1.0);

        // philipc: force weights to 1.0: required for block coupled solver
        // also arguably for stable for segregated

//         label pI = 0;
//         for (label i=0; i<allPoints.size(); i++)
//         {
//             scalar curR =  mag(allPoints[i] - p[pointI]);
// //             W[pI++] = 1.0/(curR + VSMALL);
//             W[pI++] = 1.0/(sqr(curR) + VSMALL);
//         }
//         for (label i=0; i<allMirrorPoints.size(); i++)
//         {
//             scalar curR =  mag(allMirrorPoints[i] - p[pointI]);
// //             W[pI++] = 1.0/(curR + VSMALL);
//             W[pI++] = 1.0/(sqr(curR) + VSMALL);
//         }

        weights.set(pointI, new scalarField(W));
    }
}


void newLeastSquaresVolPointInterpolation::makeOrigins() const
{
    if (debug)
    {
        Info<< "newLeastSquaresVolPointInterpolation::makeOrigin() : "
            << "constructing local origins"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!originsPtr_.empty())
    {
        FatalErrorIn
        (
            "newLeastSquaresVolPointInterpolation::makeOrigin() const"
        )
            << "local origins already exist"
            << abort(FatalError);
    }

    originsPtr_.set(new vectorField(mesh().points().size(), vector::zero));
    vectorField& origins = originsPtr_();

    refLPtr_.set(new scalarField(mesh().points().size(), 0));
    scalarField& refL = refLPtr_();

    const FieldField<Field, scalar>& w = weights();

    const vectorField& p = mesh().points();
    const vectorField& C = mesh().cellCentres();
    const vectorField& Cf = mesh().faceCentres();

    const labelListList& ptCells = mesh().pointCells();
    const labelListList& ptBndFaces = pointBndFaces();
    const labelListList& ptCyclicFaces = pointCyclicFaces();
    const labelListList& ptCyclicBndFaces = pointCyclicBndFaces();
    const labelListList& ptProcFaces = pointProcFaces();

//     const Map<vectorField> ptNgbProcBndFaceCentres =
//         pointNgbProcBndFaceCentres();

    const Map<vectorField> gPtNgbProcBndFaceCentres =
        globalPointNgbProcBndFaceCentres();
    const Map<vectorField> gPtNgbProcCellCentres =
        globalPointNgbProcCellCentres();

    const Map<List<labelPair> >& ptProcCells = pointProcCells();

    const List<List<labelPair> >& ptProcBndFaces = pointProcBndFaces();
//     const Map<List<labelPair> >& ptProcBndFaces = pointProcBndFaces();

    const FieldField<Field, vector>& procCentres = procCellCentres();
    const FieldField<Field, vector>& procBndFaceCent = procBndFaceCentres();

    forAll(origins, pointI)
    {
        const labelList& interpCells = ptCells[pointI];
        const labelList& interpBndFaces = ptBndFaces[pointI];
        const labelList& interpCyclicFaces = ptCyclicFaces[pointI];
        const labelList& interpCyclicBndFaces = ptCyclicBndFaces[pointI];
        const labelList& interpProcFaces = ptProcFaces[pointI];

//         vectorField interpNgbProcBndFaceCentres(0);

//         // Boundary faces from neighbour processors
//         if (ptNgbProcBndFaceCentres.found(pointI))
//         {
//              interpNgbProcBndFaceCentres =
//                 ptNgbProcBndFaceCentres[pointI];
//         }

        vectorField glInterpNgbProcBndFaceCentres(0);

        // Boundary faces from neighbour processors
        if (gPtNgbProcBndFaceCentres.found(pointI))
        {
             glInterpNgbProcBndFaceCentres =
                gPtNgbProcBndFaceCentres[pointI];
        }

        vectorField glInterpNgbProcCellCentres(0);

        // Boundary faces from neighbour processors
        if (gPtNgbProcCellCentres.found(pointI))
        {
             glInterpNgbProcCellCentres =
                gPtNgbProcCellCentres[pointI];
        }

        vectorField interpNgbProcCellCentres(0);

        if (ptProcCells.found(pointI))
        {
            const List<labelPair>& pc = ptProcCells[pointI];

            interpNgbProcCellCentres.setSize(pc.size());

            forAll(pc, cI)
            {
                interpNgbProcCellCentres[cI] =
                    procCentres
                    [
                        pc[cI].first()
                    ]
                    [
                        pc[cI].second()
                    ];
            }
        }

        vectorField interpNgbProcBndFaceCentres(0);

        if (ptProcBndFaces[pointI].size())
//         if (ptProcBndFaces.found(pointI))
        {
            const List<labelPair>& pf = ptProcBndFaces[pointI];

            interpNgbProcBndFaceCentres.setSize(pf.size());

            forAll(pf, fI)
            {
                interpNgbProcBndFaceCentres[fI] =
                    procBndFaceCent
                    [
                        pf[fI].first()
                    ]
                    [
                        pf[fI].second()
                    ];
            }
        }

        vectorField allPoints
        (
            interpCells.size()
          + interpBndFaces.size()
          + interpCyclicFaces.size()
          + interpCyclicBndFaces.size()
          + interpProcFaces.size()
//           + interpNgbProcBndFaceCentres.size()
          + glInterpNgbProcBndFaceCentres.size()
          + glInterpNgbProcCellCentres.size()
          + interpNgbProcCellCentres.size()
          + interpNgbProcBndFaceCentres.size(),
            vector::zero
        );

        label pointID = 0;

        // Cells
        for (label i=0; i<interpCells.size(); i++)
        {
            allPoints[pointID++] = C[interpCells[i]];
        }

        // Boundary faces
        for (label i=0; i<interpBndFaces.size(); i++)
        {
            label faceID = interpBndFaces[i];

            allPoints[pointID++] = Cf[faceID];
        }

        // Cyclic boundary faces
        for (label i=0; i<interpCyclicFaces.size(); i++)
        {
            label faceID = interpCyclicFaces[i];
            label patchID = mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            const unallocLabelList& faceCells =
                mesh().boundary()[patchID].faceCells();

            const cyclicFvPatch& cycPatch =
                refCast<const cyclicFvPatch>(mesh().boundary()[patchID]);

            label sizeby2 = faceCells.size()/2;

            if (cycPatch.parallel())
            {
                if (localFaceID < sizeby2)
                {
                    vector deltaNgb = //dni
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID + sizeby2
                        ]
                      - C[faceCells[localFaceID + sizeby2]];

                    vector deltaOwn = //ddi
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID
                        ]
                      - C[faceCells[localFaceID]];

                    vector delta = deltaOwn - deltaNgb;

                    allPoints[pointID++] = C[faceCells[localFaceID]] + delta;
                }
                else
                {
                    vector deltaNgb = //ddi
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID - sizeby2
                        ]
                      - C[faceCells[localFaceID - sizeby2]];

                    vector deltaOwn = //dni
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID
                        ]
                      - C[faceCells[localFaceID]];

                    vector delta = deltaNgb - deltaOwn;
                    delta *= -1;

                    allPoints[pointID++] = C[faceCells[localFaceID]] + delta;
                }
            }
            else
            {
                if (localFaceID < sizeby2)
                {
                    vector deltaNgb = // dni
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID + sizeby2
                        ]
                      - C[faceCells[localFaceID + sizeby2]];

                    vector deltaOwn = //ddi
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID
                        ]
                      - C[faceCells[localFaceID]];

                    vector delta =
                        deltaOwn
                      - transform(cycPatch.forwardT()[0], deltaNgb);

                    allPoints[pointID++] = C[faceCells[localFaceID]] + delta;

//                     Info << "own: " <<
//                         mesh().Cf().boundaryField()[patchID]
//                         [
//                             localFaceID
//                         ] << ", "
//                         << C[faceCells[localFaceID]]
//                         << ", " << C[faceCells[localFaceID]] + delta
//                         << endl;
                }
                else
                {
                    vector deltaNgb = //ddi
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID - sizeby2
                        ]
                      - C[faceCells[localFaceID - sizeby2]];

                    vector deltaOwn = //dni
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID
                        ]
                      - C[faceCells[localFaceID]];

                    vector delta =
                        deltaNgb
                      - transform(cycPatch.forwardT()[0], deltaOwn);

                    delta = -transform(cycPatch.reverseT()[0], delta);

                    allPoints[pointID++] = C[faceCells[localFaceID]] + delta;

//                     Info << "ngb: " <<
//                         mesh().Cf().boundaryField()[patchID]
//                         [
//                             localFaceID
//                         ] << ", "
//                         << C[faceCells[localFaceID]]
//                         << ", " << C[faceCells[localFaceID]] + delta
//                         << endl;
                }
            }

//             label faceID = interpCyclicFaces[i];
//             label patchID = mesh().boundaryMesh().whichPatch(faceID);

//             label start = mesh().boundaryMesh()[patchID].start();
//             label localFaceID = faceID - start;

//             const unallocLabelList& faceCells =
//                 mesh().boundary()[patchID].faceCells();

//             label sizeby2 = faceCells.size()/2;

//             if (localFaceID < sizeby2)
//             {
//                 vector delta =
//                     C[faceCells[localFaceID + sizeby2]]
//                   - mesh().Cf().boundaryField()[patchID]
//                     [
//                         localFaceID + sizeby2
//                     ];

//                 allPoints[pointID++] = Cf[faceID] + delta;
//             }
//             else
//             {
//                 vector delta =
//                     C[faceCells[localFaceID - sizeby2]]
//                   - mesh().Cf().boundaryField()[patchID]
//                     [
//                         localFaceID - sizeby2
//                     ];

//                 allPoints[pointID++] = Cf[faceID] + delta;
//             }
        }

        for (label i=0; i<interpCyclicBndFaces.size(); i++)
        {
            label cycFaceID = interpCyclicFaces[0];
            label cycPatchID = mesh().boundaryMesh().whichPatch(cycFaceID);
            label cycStart = mesh().boundaryMesh()[cycPatchID].start();
            label cycLocalFaceID = cycFaceID - cycStart;

            const cyclicFvPatch& cycPatch =
                refCast<const cyclicFvPatch>(mesh().boundary()[cycPatchID]);

            const cyclicPolyPatch& cycPolyPatch =
                refCast<const cyclicPolyPatch>
                (
                    mesh().boundaryMesh()[cycPatchID]
                );

            label faceID = interpCyclicBndFaces[i];

            const edgeList& coupledPoints = cycPolyPatch.coupledPoints();
            const labelList& meshPoints = cycPolyPatch.meshPoints();

            label sizeby2 = cycPolyPatch.size()/2;

            label ngbCycPointI = -1;
            forAll(coupledPoints, pI)
            {
                if (cycLocalFaceID < sizeby2)
                {
                    if ( pointI == meshPoints[coupledPoints[pI][0]] )
                    {
                        ngbCycPointI = meshPoints[coupledPoints[pI][1]];
                        break;
                    }
                }
                else
                {
                    if ( pointI == meshPoints[coupledPoints[pI][1]] )
                    {
                        ngbCycPointI = meshPoints[coupledPoints[pI][0]];
                        break;
                    }
                }
            }

            if (cycPatch.parallel())
            {
                vector deltaNgb = //dni
                    Cf[faceID] - p[ngbCycPointI];

                allPoints[pointID++] = p[pointI] + deltaNgb;
            }
            else
            {
                vector deltaNgb = //dni
                    Cf[faceID] - p[ngbCycPointI];


                if (cycLocalFaceID < sizeby2)
                {
                    allPoints[pointID++] =
                        p[pointI]
                      + transform(cycPatch.forwardT()[0], deltaNgb);
                }
                else
                {
                    allPoints[pointID++] =
                        p[pointI]
                      + transform(cycPatch.reverseT()[0], deltaNgb);
                }
            }
        }

        // Processor boundary faces
        for (label i=0; i<interpProcFaces.size(); i++)
        {
            label faceID = interpProcFaces[i];
            label patchID = mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            allPoints[pointID++] =
                mesh().C().boundaryField()[patchID][localFaceID];
        }

//         // Boundary faces from neighbour processors
//         for (label i=0; i<interpNgbProcBndFaceCentres.size(); i++)
//         {
//             allPoints[pointID++] =
//                 interpNgbProcBndFaceCentres[i];
//         }

        // Global point bnd faces from neighbour processors
        for (label i=0; i<glInterpNgbProcBndFaceCentres.size(); i++)
        {
            allPoints[pointID++] =
                glInterpNgbProcBndFaceCentres[i];
        }

        // Global point cells from neighbour processors
        for (label i=0; i<glInterpNgbProcCellCentres.size(); i++)
        {
            allPoints[pointID++] =
                glInterpNgbProcCellCentres[i];
        }

        // Cells from neighbour processors
        for (label i=0; i<interpNgbProcCellCentres.size(); i++)
        {
            allPoints[pointID++] =
                interpNgbProcCellCentres[i];
        }

        // Boundary faces from neighbour processors
        for (label i=0; i<interpNgbProcBndFaceCentres.size(); i++)
        {
            allPoints[pointID++] =
                interpNgbProcBndFaceCentres[i];
        }

        vectorField allMirrorPoints(0);
        if (mag(mirrorPlaneTransformation()[pointI].first())>SMALL)
//         if (mirrorPlaneTransformation().found(pointI))
        {
            const vector& n = mirrorPlaneTransformation()[pointI].first();

            allMirrorPoints.setSize(allPoints.size());

            forAll(allPoints, pI)
            {
                allMirrorPoints[pI] =
                    p[pointI] + transform(I-2*n*n, (allPoints[pI]-p[pointI]));
            }
        }


        const scalarField& W = w[pointI];

        label pI = 0;
        for (label i=0; i<allPoints.size(); i++)
        {
            origins[pointI] += sqr(W[pI++])*allPoints[i];
        }
        for (label i=0; i<allMirrorPoints.size(); i++)
        {
            origins[pointI] += sqr(W[pI++])*allMirrorPoints[i];
        }

        origins[pointI] /= sum(sqr(W));

        // PC, for now, disable as it negatively affects the coupled solver
        //boundBox bb(allPoints, false);
        //refL[pointI] = mag(bb.max() - bb.min())/2;
        refL[pointI] = 1.0;
    }
}


void newLeastSquaresVolPointInterpolation::makeInvLsMatrices() const
{
    if (debug)
    {
        Info<< "newLeastSquaresVolPointInterpolation::makeInvLsMatrices() : "
            << "making least squares linear interpolation matrices"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (invLsMatrices_.size() != 0)
    {
        FatalErrorIn
        (
            "newLeastSquaresVolPointInterpolation::makeInvLsMatrices()"
        )   << "least square linear inerpolation matrices already exist"
            << abort(FatalError);
    }

    invLsMatrices_.setSize(mesh().points().size());

    const vectorField& p = mesh().points();
    const vectorField& C = mesh().cellCentres();
    const vectorField& Cf = mesh().faceCentres();

    const labelListList& ptCells = mesh().pointCells();
    const labelListList& ptBndFaces = pointBndFaces();
    const labelListList& ptCyclicFaces = pointCyclicFaces();
    const labelListList& ptCyclicBndFaces = pointCyclicBndFaces();
    const labelListList& ptProcFaces = pointProcFaces();

//     const Map<vectorField>& ptNgbProcBndFaceCentres =
//         pointNgbProcBndFaceCentres();

    const Map<vectorField>& gPtNgbProcBndFaceCentres =
        globalPointNgbProcBndFaceCentres();

    const Map<vectorField>& gPtNgbProcCellCentres =
        globalPointNgbProcCellCentres();

    const Map<List<labelPair> >& ptProcCells = pointProcCells();
    const FieldField<Field, vector>& procCentres = procCellCentres();

    const List<List<labelPair> >& ptProcBndFaces = pointProcBndFaces();
//     const Map<List<labelPair> >& ptProcBndFaces = pointProcBndFaces();
    const FieldField<Field, vector>& procBndFaceCent = procBndFaceCentres();

    const FieldField<Field, scalar>& w = weights();
    const vectorField& o = origins();

    const scalarField& L = refL();

    label nCoeffs = 3;

    scalarField D(invLsMatrices_.size(), 0);

    forAll(invLsMatrices_, pointI)
    {
        const labelList& interpCells = ptCells[pointI];
        const labelList& interpBndFaces = ptBndFaces[pointI];
        const labelList& interpCyclicFaces = ptCyclicFaces[pointI];
        const labelList& interpCyclicBndFaces = ptCyclicBndFaces[pointI];
        const labelList& interpProcFaces = ptProcFaces[pointI];

//         vectorField interpNgbProcBndFaceCentres(0);// = vectorField::null();

//         // Boundary faces from neighbour processors
//         if (ptNgbProcBndFaceCentres.found(pointI))
//         {
//              interpNgbProcBndFaceCentres =
//                 ptNgbProcBndFaceCentres[pointI];
//         }

        vectorField glInterpNgbProcBndFaceCentres(0);

        // Boundar faces from neighbour processors
        if (gPtNgbProcBndFaceCentres.found(pointI))
        {
             glInterpNgbProcBndFaceCentres =
                gPtNgbProcBndFaceCentres[pointI];
        }

        vectorField glInterpNgbProcCellCentres(0);

        // Cells from neighbour processors
        if (gPtNgbProcCellCentres.found(pointI))
        {
             glInterpNgbProcCellCentres =
                gPtNgbProcCellCentres[pointI];
        }

        vectorField interpNgbProcCellCentres(0);

        if (ptProcCells.found(pointI))
        {
            const List<labelPair>& pc = ptProcCells[pointI];

            interpNgbProcCellCentres.setSize(pc.size());

            forAll(pc, cI)
            {
                interpNgbProcCellCentres[cI] =
                    procCentres
                    [
                        pc[cI].first()
                    ]
                    [
                        pc[cI].second()
                    ];
            }
        }

        vectorField interpNgbProcBndFaceCentres(0);

        if (ptProcBndFaces[pointI].size())
//         if (ptProcBndFaces.found(pointI))
        {
            const List<labelPair>& pf = ptProcBndFaces[pointI];

            interpNgbProcBndFaceCentres.setSize(pf.size());

            forAll(pf, fI)
            {
                interpNgbProcBndFaceCentres[fI] =
                    procBndFaceCent
                    [
                        pf[fI].first()
                    ]
                    [
                        pf[fI].second()
                    ];
            }
        }

        vectorField allPoints
        (
            interpCells.size()
          + interpBndFaces.size()
          + interpCyclicFaces.size()
          + interpCyclicBndFaces.size()
          + interpProcFaces.size()
//           + interpNgbProcBndFaceCentres.size()
          + glInterpNgbProcBndFaceCentres.size()
          + glInterpNgbProcCellCentres.size()
          + interpNgbProcCellCentres.size()
          + interpNgbProcBndFaceCentres.size(),
            vector::zero
        );

        label pointID = 0;

        // Cells
        for (label i=0; i<interpCells.size(); i++)
        {
            allPoints[pointID++] = C[interpCells[i]];
        }

        // Boundary faces
        for (label i=0; i<interpBndFaces.size(); i++)
        {
            label faceID = interpBndFaces[i];

            allPoints[pointID++] = Cf[faceID];
        }

        // Cyclic boundary faces
        for (label i=0; i<interpCyclicFaces.size(); i++)
        {
            label faceID = interpCyclicFaces[i];
            label patchID = mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            const unallocLabelList& faceCells =
                mesh().boundary()[patchID].faceCells();

            const cyclicFvPatch& cycPatch =
                refCast<const cyclicFvPatch>(mesh().boundary()[patchID]);

            label sizeby2 = faceCells.size()/2;

            if (cycPatch.parallel())
            {
                if (localFaceID < sizeby2)
                {
                    vector deltaNgb = //dni
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID + sizeby2
                        ]
                      - C[faceCells[localFaceID + sizeby2]];

                    vector deltaOwn = //ddi
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID
                        ]
                      - C[faceCells[localFaceID]];

                    vector delta = deltaOwn - deltaNgb;

                    allPoints[pointID++] = C[faceCells[localFaceID]] + delta;
                }
                else
                {
                    vector deltaNgb = //ddi
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID - sizeby2
                        ]
                      - C[faceCells[localFaceID - sizeby2]];

                    vector deltaOwn = //dni
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID
                        ]
                      - C[faceCells[localFaceID]];

                    vector delta = deltaNgb - deltaOwn;
                    delta *= -1;

                    allPoints[pointID++] = C[faceCells[localFaceID]] + delta;
                }
            }
            else
            {
                if (localFaceID < sizeby2)
                {
                    vector deltaNgb = // dni
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID + sizeby2
                        ]
                      - C[faceCells[localFaceID + sizeby2]];

                    vector deltaOwn = //ddi
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID
                        ]
                      - C[faceCells[localFaceID]];

                    vector delta =
                        deltaOwn
                      - transform(cycPatch.forwardT()[0], deltaNgb);

                    allPoints[pointID++] = C[faceCells[localFaceID]] + delta;

//                     Info << "own: " <<
//                         mesh().Cf().boundaryField()[patchID]
//                         [
//                             localFaceID
//                         ] << ", "
//                         << C[faceCells[localFaceID]]
//                         << ", " << C[faceCells[localFaceID]] + delta
//                         << endl;
                }
                else
                {
                    vector deltaNgb = //ddi
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID - sizeby2
                        ]
                      - C[faceCells[localFaceID - sizeby2]];

                    vector deltaOwn = //dni
                        mesh().Cf().boundaryField()[patchID]
                        [
                            localFaceID
                        ]
                      - C[faceCells[localFaceID]];

                    vector delta =
                        deltaNgb
                      - transform(cycPatch.forwardT()[0], deltaOwn);

                    delta = -transform(cycPatch.reverseT()[0], delta);

                    allPoints[pointID++] = C[faceCells[localFaceID]] + delta;

//                     Info << "ngb: " <<
//                         mesh().Cf().boundaryField()[patchID]
//                         [
//                             localFaceID
//                         ] << ", "
//                         << C[faceCells[localFaceID]]
//                         << ", " << C[faceCells[localFaceID]] + delta
//                         << endl;
                }
            }

//             label faceID = interpCyclicFaces[i];
//             label patchID = mesh().boundaryMesh().whichPatch(faceID);

//             label start = mesh().boundaryMesh()[patchID].start();
//             label localFaceID = faceID - start;

//             const unallocLabelList& faceCells =
//                 mesh().boundary()[patchID].faceCells();

//             label sizeby2 = faceCells.size()/2;

//             if (localFaceID < sizeby2)
//             {
//                 vector delta =
//                     C[faceCells[localFaceID + sizeby2]]
//                   - mesh().Cf().boundaryField()[patchID]
//                     [
//                         localFaceID + sizeby2
//                     ];

//                 allPoints[pointID++] = Cf[faceID] + delta;
//             }
//             else
//             {
//                 vector delta =
//                     C[faceCells[localFaceID - sizeby2]]
//                   - mesh().Cf().boundaryField()[patchID]
//                     [
//                         localFaceID - sizeby2
//                     ];

//                 allPoints[pointID++] = Cf[faceID] + delta;
//             }
        }

        for (label i=0; i<interpCyclicBndFaces.size(); i++)
        {
            label cycFaceID = interpCyclicFaces[0];
            label cycPatchID = mesh().boundaryMesh().whichPatch(cycFaceID);
            label cycStart = mesh().boundaryMesh()[cycPatchID].start();
            label cycLocalFaceID = cycFaceID - cycStart;

            const cyclicFvPatch& cycPatch =
                refCast<const cyclicFvPatch>(mesh().boundary()[cycPatchID]);

            const cyclicPolyPatch& cycPolyPatch =
                refCast<const cyclicPolyPatch>
                (
                    mesh().boundaryMesh()[cycPatchID]
                );

            label faceID = interpCyclicBndFaces[i];

            const edgeList& coupledPoints = cycPolyPatch.coupledPoints();
            const labelList& meshPoints = cycPolyPatch.meshPoints();

            label sizeby2 = cycPolyPatch.size()/2;

            label ngbCycPointI = -1;
            forAll(coupledPoints, pI)
            {
                if (cycLocalFaceID < sizeby2)
                {
                    if ( pointI == meshPoints[coupledPoints[pI][0]] )
                    {
                        ngbCycPointI = meshPoints[coupledPoints[pI][1]];
                        break;
                    }
                }
                else
                {
                    if ( pointI == meshPoints[coupledPoints[pI][1]] )
                    {
                        ngbCycPointI = meshPoints[coupledPoints[pI][0]];
                        break;
                    }
                }
            }

            if (cycPatch.parallel())
            {
                vector deltaNgb = //dni
                    Cf[faceID] - p[ngbCycPointI];

                allPoints[pointID++] = p[pointI] + deltaNgb;
            }
            else
            {
                vector deltaNgb = //dni
                    Cf[faceID] - p[ngbCycPointI];


                if (cycLocalFaceID < sizeby2)
                {
                    allPoints[pointID++] =
                        p[pointI]
                      + transform(cycPatch.transform().R()[0], deltaNgb);
                }
                else
                {
                    allPoints[pointID++] =
                        p[pointI]
                      + transform(cycPatch.transform().R()[0].T(), deltaNgb);
                }
            }
        }

        // Processor boundary faces
        for (label i=0; i<interpProcFaces.size(); i++)
        {
            label faceID = interpProcFaces[i];
            label patchID = mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            allPoints[pointID++] =
                mesh().C().boundaryField()[patchID][localFaceID];
        }

//         // Boundary faces from neighbour processors
//         for (label i=0; i<interpNgbProcBndFaceCentres.size(); i++)
//         {
//             allPoints[pointID++] =
//                 interpNgbProcBndFaceCentres[i];
//         }

        // Global point bnd faces from neighbour processors
        for (label i=0; i<glInterpNgbProcBndFaceCentres.size(); i++)
        {
            allPoints[pointID++] =
                glInterpNgbProcBndFaceCentres[i];
        }

        // Global point cells from neighbour processors
        for (label i=0; i<glInterpNgbProcCellCentres.size(); i++)
        {
            allPoints[pointID++] =
                glInterpNgbProcCellCentres[i];
        }

        // Cells from neighbour processors
        for (label i=0; i<interpNgbProcCellCentres.size(); i++)
        {
            allPoints[pointID++] =
                interpNgbProcCellCentres[i];
        }

        // Bnd faces from neighbour processors
        for (label i=0; i<interpNgbProcBndFaceCentres.size(); i++)
        {
            allPoints[pointID++] =
                interpNgbProcBndFaceCentres[i];
        }

        vectorField allMirrorPoints(0);
        if (mag(mirrorPlaneTransformation()[pointI].first()) > SMALL)
//         if (mirrorPlaneTransformation().found(pointI))
        {
            const vector& n = mirrorPlaneTransformation()[pointI].first();

            allMirrorPoints.setSize(allPoints.size());

            forAll(allPoints, pI)
            {
                allMirrorPoints[pI] =
                    p[pointI] + transform(I-2*n*n, (allPoints[pI]-p[pointI]));
            }
        }

        if (allPoints.size() + allMirrorPoints.size() < nCoeffs)
        {
            Pout << pointI << ", "
                << interpCells.size() << ", "
                << interpBndFaces.size() << ", "
                << interpCyclicFaces.size() << ", "
                << interpCyclicBndFaces.size() << ", "
                << interpProcFaces.size() << ", "
                << interpNgbProcCellCentres.size() << endl;

            FatalErrorIn
            (
                "newLeastSquaresVolPointInterpolation::makeInvInvMatrices()"
            )   << "allPoints.size() < " << nCoeffs << " : "
                << allPoints.size() + allMirrorPoints.size()
                << abort(FatalError);
        }

        // Weights
        // scalarField W(allPoints.size(), 1.0);
        // scalar sumW = 0;
        // for (label i=0; i<allPoints.size(); i++)
        // {
        //     scalar curR =  mag(allPoints[i] - p[pointI]);
        //     W[i] = 1.0/(sqr(curR) + VSMALL);
        //     sumW += W[i];
        // }
        // W /= sumW;

        const scalarField& W = w[pointI];

        invLsMatrices_.set
        (
            pointI,
            new scalarRectangularMatrix
            (
                nCoeffs,
                allPoints.size() + allMirrorPoints.size(),
                0.0
            )
        );
        scalarRectangularMatrix& curMatrix = invLsMatrices_[pointI];

        scalarRectangularMatrix M
        (
            allPoints.size() + allMirrorPoints.size(),
            nCoeffs,
            0.0
        );

        label pI = 0;
        for (label i=0; i<allPoints.size(); i++)
        {
            scalar X = (allPoints[i].x() - o[pointI].x())/L[pointI];
            scalar Y = (allPoints[i].y() - o[pointI].y())/L[pointI];
            scalar Z = (allPoints[i].z() - o[pointI].z())/L[pointI];

            M[pI][0] = X;
            M[pI][1] = Y;
            M[pI][2] = Z;
            pI++;
        }
        for (label i=0; i<allMirrorPoints.size(); i++)
        {
            scalar X = (allMirrorPoints[i].x() - o[pointI].x())/L[pointI];
            scalar Y = (allMirrorPoints[i].y() - o[pointI].y())/L[pointI];
            scalar Z = (allMirrorPoints[i].z() - o[pointI].z())/L[pointI];

            M[pI][0] = X;
            M[pI][1] = Y;
            M[pI][2] = Z;
            pI++;
        }

        // Note: the definition of n() and m() are flipped in foam-extend-4.0
        // and OpenFOAM... wtf
        // Applying weights
        for (label i=0; i<M.m(); i++)
        {
            for (label j=0; j<M.n(); j++)
            {
                M[i][j] *= W[i];
            }
        }

        //         SVD svd(M, SMALL);

//         for (label i=0; i<svd.VSinvUt().n(); i++)
//         {
//             for (label j=0; j<svd.VSinvUt().m(); j++)
//             {
//                 curMatrix[i][j] = svd.VSinvUt()[i][j]*W[j];
//             }
//         }

//         scalarSquareMatrix lsM(nCoeffs, 0.0);

        tensor lsM = tensor::zero;

        for (label i=0; i<3; i++)
        {
            for (label j=0; j<3; j++)
            {
                for (label k=0; k<M.m(); k++)
                {
//                     lsM[i][j] += M[k][i]*M[k][j];
                    lsM(i,j) += M[k][i]*M[k][j];
                }
            }
        }

        // Calculate matrix norm
        scalar maxRowSum = 0.0;
        for (label i=0; i<3; i++)
        {
            scalar curRowSum = 0.0;

            for (label j=0; j<3; j++)
            {
//                 curRowSum += lsM[i][j];
                curRowSum += lsM(i,j);
            }
            if(curRowSum > maxRowSum)
            {
                maxRowSum = curRowSum;
            }
        }

        // Calculate inverse
        D[pointI] = det(lsM);

        // PC, for now, disable as it negatively affects the coupled solver
        // if (mag(D[pointI]) > SMALL)
        {
            tensor invLsM = hinv(lsM);
            //tensor invLsM = inv(lsM);

            for (label i=0; i<3; i++)
            {
                for (label j=0; j<M.m(); j++)
                {
                    for (label k=0; k<3; k++)
                    {
                        curMatrix[i][j] += invLsM(i,k)*M[j][k]*W[j];
                    }
                }
            }
        }
        // else
        // {
        //     Pout<< "Det: " << D[pointI] << endl;

        //     for (label i=0; i<3; i++)
        //     {
        //         for (label j=0; j<M.n(); j++)
        //         {
        //             curMatrix[i][j] = 0;
        //         }
        //     }
        // }
    }
}


void newLeastSquaresVolPointInterpolation::
makeMirrorPlaneTransformation() const
{
    if (debug)
    {
        Info<< "newLeastSquaresVolPointInterpolation::"
            << "makeMirrorPlaneTransformation() : "
            << "constructing mirror plane normals and transformation tensors"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!mirrorPlaneTransformationPtr_.empty())
    {
        FatalErrorIn
        (
            "newLeastSquaresVolPointInterpolation::"
            "makeMirrorPlaneTransformation() const"
        )
            << "Mirror plane normals and transformation tensors already exist"
                << abort(FatalError);
    }

    mirrorPlaneTransformationPtr_.set
    (
        new List<Tuple2<vector, tensor> >
        (
            mesh().points().size(),
            Tuple2<vector, tensor>(vector::zero, tensor::zero)
        )
    );
    List<Tuple2<vector, tensor> >& mirrorPlaneTransformation =
        mirrorPlaneTransformationPtr_();

//     mirrorPlaneTransformationPtr_ = new Map<Tuple2<vector, tensor> >();
//     Map<Tuple2<vector, tensor> >& mirrorPlaneTransformation =
//         *mirrorPlaneTransformationPtr_;


    forAll(mesh().boundaryMesh(), patchI)
    {
        if
        (
            (
                mesh().boundaryMesh()[patchI].type()
             == emptyPolyPatch::typeName
            )
        )
        {
            const labelList& meshPoints =
                mesh().boundaryMesh()[patchI].meshPoints();

            const vectorField& pointNormals =
                mesh().boundaryMesh()[patchI].pointNormals();

            forAll(meshPoints, pointI)
            {
                mirrorPlaneTransformation[meshPoints[pointI]] =
                    Tuple2<vector, tensor>
                    (
                        pointNormals[pointI],
                        I
                    );

//                 mirrorPlaneTransformation.insert
//                 (
//                     meshPoints[pointI],
//                     Tuple2<vector, tensor>
//                     (
//                         pointNormals[pointI],
//                         I
//                     )
//                 );
            }
        }
        else if
        (
            mesh().boundaryMesh()[patchI].type()
         == wedgePolyPatch::typeName
        )
        {
            const labelList& meshPoints =
                mesh().boundaryMesh()[patchI].meshPoints();

            const vectorField& pointNormals =
                mesh().boundaryMesh()[patchI].pointNormals();

            const wedgePolyPatch& wedge =
                refCast<const wedgePolyPatch>(mesh().boundaryMesh()[patchI]);

            forAll(meshPoints, pointI)
            {
                if (!pointAxisEdges().found(meshPoints[pointI]))
                {
                    mirrorPlaneTransformation[meshPoints[pointI]] =
                        Tuple2<vector, tensor>
                        (
                            pointNormals[pointI],
                            wedge.cellT()
                        );
                }
//                 mirrorPlaneTransformation.insert
//                 (
//                     meshPoints[pointI],
//                     Tuple2<vector, tensor>
//                     (
//                         pointNormals[pointI],
//                         wedge.cellT()
//                     )
//                 );
            }
        }
    }
}

tensor newLeastSquaresVolPointInterpolation::hinv(const tensor& t) const
{
    static const scalar hinvLarge = 1e10;
    static const scalar hinvSmall = 1e-10;

    if (det(t) > hinvSmall)
    {
        return inv(t);
    }
    else
    {
        vector eig = eigenValues(t);
        tensor eigVecs = eigenVectors(t);

        tensor zeroInv = tensor::zero;

        // Test if all eigen values are zero.
        // If this happens then eig.z() = SMALL, and hinv(t)
        // returns a zero tensor.
        // Jovani Favero, 18/Nov/2009
        // Further bug fix: replace > with == and add SMALL to zeroInv
        // Dominik Christ, 7/Aug/2012
        if (mag(eig.z()) == hinvLarge*mag(eig.z()))
        {
            // Three zero eigen values (z is largest in magnitude).
            // Return zero inverse
            return zeroInv;
        }

        // Compare smaller eigen values and if out of range of large
        // consider them singular

        if (mag(eig.z()) > hinvLarge*mag(eig.x()))
        {
            // Make a tensor out of symmTensor sqr.  HJ, 24/Oct/2009
            zeroInv += tensor(sqr(eigVecs.x()));
        }

        if (mag(eig.z()) > hinvLarge*mag(eig.y()))
        {
            // Make a tensor out of symmTensor sqr.  HJ, 24/Oct/2009
            zeroInv += tensor(sqr(eigVecs.y()));
        }

        return inv(t + zeroInv) - zeroInv;
    }
}

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

newLeastSquaresVolPointInterpolation::newLeastSquaresVolPointInterpolation
(
    const fvMesh& vm
)
:
    MeshObject<fvMesh, Foam::UpdateableMeshObject, newLeastSquaresVolPointInterpolation>(vm),
    pointBndFacesPtr_(),
    pointCyclicFacesPtr_(),
    pointCyclicBndFacesPtr_(),
    pointCyclicGgiFacesPtr_(),
    pointCyclicGgiBndFacesPtr_(),
    pointProcFacesPtr_(),
    axisEdgesPtr_(),
    pointAxisEdgesPtr_(),
    globalPointNgbProcBndFaceCentresPtr_(),
    globalPointNgbProcCellCentresPtr_(),
    procBndFacesPtr_(),
    procBndFaceCentresPtr_(),
    pointProcBndFacesPtr_(),
    procCellsPtr_(),
    pointProcCellsPtr_(),
    procCellCentresPtr_(),
    weightsPtr_(),
    originsPtr_(),
    mirrorPlaneTransformationPtr_(),
    invLsMatrices_(0),
    refLPtr_(),
    processorBoundariesExist_(false)
{
    if (Pstream::parRun())
    {
        // Check there are no processor boundaries
        forAll(mesh().boundaryMesh(), patchI)
        {
            if (mesh().boundaryMesh()[patchI].type() == "processor")
            {
                processorBoundariesExist_ = true;
                break;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

newLeastSquaresVolPointInterpolation::~newLeastSquaresVolPointInterpolation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool newLeastSquaresVolPointInterpolation::movePoints()
{
    // Clear all fields that depend on the mesh coordinates
    globalPointNgbProcBndFaceCentresPtr_.clear();
    globalPointNgbProcCellCentresPtr_.clear();
    procBndFaceCentresPtr_.clear();
    procCellCentresPtr_.clear();
    weightsPtr_.clear();
    originsPtr_.clear();
    mirrorPlaneTransformationPtr_.clear();
    invLsMatrices_.clear();
    refLPtr_.clear();

    return true;
}


void newLeastSquaresVolPointInterpolation::updateMesh(const mapPolyMesh&)
{
    // Clear all fields that depend on the mesh addressing
    pointBndFacesPtr_.clear();
    pointCyclicFacesPtr_.clear();
    pointCyclicBndFacesPtr_.clear();
    pointCyclicGgiFacesPtr_.clear();
    pointCyclicGgiBndFacesPtr_.clear();
    pointProcFacesPtr_.clear();
    axisEdgesPtr_.clear();
    pointAxisEdgesPtr_.clear();
    globalPointNgbProcBndFaceCentresPtr_.clear();
    globalPointNgbProcCellCentresPtr_.clear();
    procBndFacesPtr_.clear();
    procBndFaceCentresPtr_.clear();
    pointProcBndFacesPtr_.clear();
    procCellsPtr_.clear();
    pointProcCellsPtr_.clear();
    procCellCentresPtr_.clear();
    weightsPtr_.clear();
    originsPtr_.clear();
    mirrorPlaneTransformationPtr_.clear();
    invLsMatrices_.clear();
    refLPtr_.clear();
}

const labelListList& newLeastSquaresVolPointInterpolation::pointBndFaces() const
{
    if (pointBndFacesPtr_.empty())
    {
        makePointFaces();
    }

    return pointBndFacesPtr_();
}

const labelListList& newLeastSquaresVolPointInterpolation
::pointCyclicFaces() const
{
    if (pointCyclicFacesPtr_.empty())
    {
        makePointFaces();
    }

    return pointCyclicFacesPtr_();
}

const labelListList& newLeastSquaresVolPointInterpolation
::pointCyclicBndFaces() const
{
    if (pointCyclicBndFacesPtr_.empty())
    {
        makePointFaces();
    }

    return pointCyclicBndFacesPtr_();
}

const labelListList& newLeastSquaresVolPointInterpolation
::pointCyclicGgiFaces() const
{
    if (pointCyclicGgiFacesPtr_.empty())
    {
        makePointFaces();
    }

    return pointCyclicGgiFacesPtr_();
}

const labelListList& newLeastSquaresVolPointInterpolation
::pointCyclicGgiBndFaces() const
{
    if (pointCyclicGgiBndFacesPtr_.empty())
    {
        makePointFaces();
    }

    return pointCyclicGgiBndFacesPtr_();
}

const labelList& newLeastSquaresVolPointInterpolation
::axisEdges() const
{
    if (axisEdgesPtr_.empty())
    {
        makeAxisEdges();
    }

    return axisEdgesPtr_();
}

const Map<labelList>&
newLeastSquaresVolPointInterpolation::pointAxisEdges() const
{
    if (pointAxisEdgesPtr_.empty())
    {
        makePointAxisEdges();
    }

    return pointAxisEdgesPtr_();
}

// const Map<Field<vector> >&
// newLeastSquaresVolPointInterpolation::pointNgbProcBndFaceCentres() const
// {
//     if (!pointNgbProcBndFaceCentresPtr_)
//     {
//         makePointNgbProcBndFaceCentres();
//     }

//     return *pointNgbProcBndFaceCentresPtr_;
// }

const Map<Field<vector> >&
newLeastSquaresVolPointInterpolation::globalPointNgbProcBndFaceCentres() const
{
    if (globalPointNgbProcBndFaceCentresPtr_.empty())
    {
        makeGlobalPointNgbProcBndFaceCentres();
    }

    return globalPointNgbProcBndFaceCentresPtr_();
}

const Map<Field<vector> >&
newLeastSquaresVolPointInterpolation::globalPointNgbProcCellCentres() const
{
    if (globalPointNgbProcCellCentresPtr_.empty())
    {
        makeGlobalPointNgbProcCellCentres();
    }

    return globalPointNgbProcCellCentresPtr_();
}

const labelListList& newLeastSquaresVolPointInterpolation::pointProcFaces() const
{
    if (pointProcFacesPtr_.empty())
    {
        makePointFaces();
    }

    return pointProcFacesPtr_();
}

const labelListList& newLeastSquaresVolPointInterpolation::procBndFaces() const
{
    if (procBndFacesPtr_.empty())
    {
        makeProcBndFaces();
    }

    return procBndFacesPtr_();
}

const FieldField<Field, vector>&
newLeastSquaresVolPointInterpolation::procBndFaceCentres() const
{
    if (procBndFaceCentresPtr_.empty())
    {
        makeProcBndFaceCentres();
    }

    return procBndFaceCentresPtr_();
}

const List<List<labelPair> >& newLeastSquaresVolPointInterpolation::
pointProcBndFaces() const
{
    if (pointProcBndFacesPtr_.empty())
    {
        makeProcBndFaces();
    }

    return pointProcBndFacesPtr_();
}

const labelListList& newLeastSquaresVolPointInterpolation::procCells() const
{
    if (procCellsPtr_.empty())
    {
        makeProcCells();
    }

    return procCellsPtr_();
}

const Map<List<labelPair> >& newLeastSquaresVolPointInterpolation::
pointProcCells() const
{
    if (pointProcCellsPtr_.empty())
    {
        makeProcCells();
    }

    return pointProcCellsPtr_();
}

const FieldField<Field, vector>&
newLeastSquaresVolPointInterpolation::procCellCentres() const
{
    if (procCellCentresPtr_.empty())
    {
        makeProcCellCentres();
    }

    return procCellCentresPtr_();
}

const FieldField<Field, scalar>&
newLeastSquaresVolPointInterpolation::weights() const
{
    if (weightsPtr_.empty())
    {
        makeWeights();
    }

    return weightsPtr_();
}

const vectorField& newLeastSquaresVolPointInterpolation::origins() const
{
    if (originsPtr_.empty())
    {
        makeOrigins();
    }

    return originsPtr_();
}

const scalarField& newLeastSquaresVolPointInterpolation::refL() const
{
    if (refLPtr_.empty())
    {
        makeOrigins();
    }

    return refLPtr_();
}

const List<Tuple2<vector, tensor> >& newLeastSquaresVolPointInterpolation::
mirrorPlaneTransformation() const
{
    if (mirrorPlaneTransformationPtr_.empty())
    {
        makeMirrorPlaneTransformation();
    }

    return mirrorPlaneTransformationPtr_();
}

const PtrList<scalarRectangularMatrix>&
newLeastSquaresVolPointInterpolation::invLsMatrices() const
{
    label size = invLsMatrices_.size();

    reduce(size, maxOp<label>());

    if (size == 0)
    {
        makeInvLsMatrices();
    }

    return invLsMatrices_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
