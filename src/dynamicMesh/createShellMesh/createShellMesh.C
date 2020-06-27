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

#include "createShellMesh.H"
#include "polyTopoChange.H"
#include "meshTools.H"
#include "mapPolyMesh.H"
#include "polyAddPoint.H"
#include "polyAddFace.H"
#include "polyModifyFace.H"
#include "polyAddCell.H"
#include "labelPair.H"
#include "indirectPrimitivePatch.H"
#include "mapDistribute.H"
#include "globalMeshData.H"
#include "PatchTools.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(createShellMesh, 0);

template<>
class minEqOp<labelPair>
{
public:
    void operator()(labelPair& x, const labelPair& y) const
    {
        x[0] = min(x[0], y[0]);
        x[1] = min(x[1], y[1]);
    }
};
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Synchronise edges
void Foam::createShellMesh::syncEdges
(
    const globalMeshData& globalData,

    const labelList& patchEdges,
    const labelList& coupledEdges,
    const PackedBoolList& sameEdgeOrientation,
    const bool syncNonCollocated,

    PackedBoolList& isChangedEdge,
    DynamicList<label>& changedEdges,
    labelPairList& allEdgeData
)
{
    const mapDistribute& map = globalData.globalEdgeSlavesMap();
    const PackedBoolList& cppOrientation = globalData.globalEdgeOrientation();

    // Convert patch-edge data into cpp-edge data
    labelPairList cppEdgeData
    (
        map.constructSize(),
        labelPair(labelMax, labelMax)
    );

    forAll(patchEdges, i)
    {
        label patchEdgeI = patchEdges[i];
        label coupledEdgeI = coupledEdges[i];

        if (isChangedEdge[patchEdgeI])
        {
            const labelPair& data = allEdgeData[patchEdgeI];

            // Patch-edge data needs to be converted into coupled-edge data
            // (optionally flipped) and consistent in orientation with
            // other coupled edge (optionally flipped)
            if (sameEdgeOrientation[i] == cppOrientation[coupledEdgeI])
            {
                cppEdgeData[coupledEdgeI] = data;
            }
            else
            {
                cppEdgeData[coupledEdgeI] = labelPair(data[1], data[0]);
            }
        }
    }

    // Synchronise
    globalData.syncData
    (
        cppEdgeData,
        globalData.globalEdgeSlaves(),
        (
            syncNonCollocated
          ? globalData.globalEdgeTransformedSlaves()
          : labelListList(globalData.globalEdgeSlaves().size())
        ),
        map,
        minEqOp<labelPair>()
    );

    // Back from cpp-edge to patch-edge data
    forAll(patchEdges, i)
    {
        label patchEdgeI = patchEdges[i];
        label coupledEdgeI = coupledEdges[i];

        if (cppEdgeData[coupledEdgeI] != labelPair(labelMax, labelMax))
        {
            const labelPair& data = cppEdgeData[coupledEdgeI];

            if (sameEdgeOrientation[i] == cppOrientation[coupledEdgeI])
            {
                allEdgeData[patchEdgeI] = data;
            }
            else
            {
                allEdgeData[patchEdgeI] = labelPair(data[1], data[0]);
            }

            if (!isChangedEdge[patchEdgeI])
            {
                changedEdges.append(patchEdgeI);
                isChangedEdge[patchEdgeI] = true;
            }
        }
    }
}


void Foam::createShellMesh::calcPointRegions
(
    const globalMeshData& globalData,
    const primitiveFacePatch& patch,
    const PackedBoolList& nonManifoldEdge,
    const bool syncNonCollocated,

    faceList& pointGlobalRegions,
    faceList& pointLocalRegions,
    labelList& localToGlobalRegion
)
{
    const indirectPrimitivePatch& cpp = globalData.coupledPatch();

    // Calculate correspondence between patch and globalData.coupledPatch.
    labelList patchEdges;
    labelList coupledEdges;
    PackedBoolList sameEdgeOrientation;
    PatchTools::matchEdges
    (
        cpp,
        patch,

        coupledEdges,
        patchEdges,
        sameEdgeOrientation
    );


    // Initial unique regions
    // ~~~~~~~~~~~~~~~~~~~~~~
    // These get merged later on across connected edges.

    // 1. Count
    label nMaxRegions = 0;
    forAll(patch.localFaces(), facei)
    {
        const face& f = patch.localFaces()[facei];
        nMaxRegions += f.size();
    }

    const globalIndex globalRegions(nMaxRegions);

    // 2. Assign unique regions
    label nRegions = 0;

    pointGlobalRegions.setSize(patch.size());
    forAll(pointGlobalRegions, facei)
    {
        const face& f = patch.localFaces()[facei];
        labelList& pRegions = pointGlobalRegions[facei];
        pRegions.setSize(f.size());
        forAll(pRegions, fp)
        {
            pRegions[fp] = globalRegions.toGlobal(nRegions++);
        }
    }


    DynamicList<label> changedEdges(patch.nEdges());
    labelPairList allEdgeData(patch.nEdges(), labelPair(labelMax, labelMax));
    PackedBoolList isChangedEdge(patch.nEdges());


    // Fill initial seed
    // ~~~~~~~~~~~~~~~~~

    forAll(patch.edgeFaces(), edgeI)
    {
        if (!nonManifoldEdge[edgeI])
        {
            // Take over value from one face only.
            const edge& e = patch.edges()[edgeI];
            label facei = patch.edgeFaces()[edgeI][0];
            const face& f = patch.localFaces()[facei];

            label fp0 = findIndex(f, e[0]);
            label fp1 = findIndex(f, e[1]);
            allEdgeData[edgeI] = labelPair
            (
                pointGlobalRegions[facei][fp0],
                pointGlobalRegions[facei][fp1]
            );
            if (!isChangedEdge[edgeI])
            {
                changedEdges.append(edgeI);
                isChangedEdge[edgeI] = true;
            }
        }
    }


    syncEdges
    (
        globalData,

        patchEdges,
        coupledEdges,
        sameEdgeOrientation,
        syncNonCollocated,

        isChangedEdge,
        changedEdges,
        allEdgeData
    );


    // Edge-Face-Edge walk across patch
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Across edge minimum regions win

    while (true)
    {
        // From edge to face
        // ~~~~~~~~~~~~~~~~~

        DynamicList<label> changedFaces(patch.size());
        PackedBoolList isChangedFace(patch.size());

        forAll(changedEdges, changedI)
        {
            label edgeI = changedEdges[changedI];
            const labelPair& edgeData = allEdgeData[edgeI];

            const edge& e = patch.edges()[edgeI];
            const labelList& eFaces = patch.edgeFaces()[edgeI];

            forAll(eFaces, i)
            {
                label facei = eFaces[i];
                const face& f = patch.localFaces()[facei];

                // Combine edgeData with face data
                label fp0 = findIndex(f, e[0]);
                if (pointGlobalRegions[facei][fp0] > edgeData[0])
                {
                    pointGlobalRegions[facei][fp0] = edgeData[0];
                    if (!isChangedFace[facei])
                    {
                        isChangedFace[facei] = true;
                        changedFaces.append(facei);
                    }
                }

                label fp1 = findIndex(f, e[1]);
                if (pointGlobalRegions[facei][fp1] > edgeData[1])
                {
                    pointGlobalRegions[facei][fp1] = edgeData[1];
                    if (!isChangedFace[facei])
                    {
                        isChangedFace[facei] = true;
                        changedFaces.append(facei);
                    }
                }
            }
        }


        label nChangedFaces = returnReduce(changedFaces.size(), sumOp<label>());
        if (nChangedFaces == 0)
        {
            break;
        }


        // From face to edge
        // ~~~~~~~~~~~~~~~~~

        isChangedEdge = false;
        changedEdges.clear();

        forAll(changedFaces, i)
        {
            label facei = changedFaces[i];
            const face& f = patch.localFaces()[facei];
            const labelList& fEdges = patch.faceEdges()[facei];

            forAll(fEdges, fp)
            {
                label edgeI = fEdges[fp];

                if (!nonManifoldEdge[edgeI])
                {
                    const edge& e = patch.edges()[edgeI];
                    label fp0 = findIndex(f, e[0]);
                    label region0 = pointGlobalRegions[facei][fp0];
                    label fp1 = findIndex(f, e[1]);
                    label region1 = pointGlobalRegions[facei][fp1];

                    if
                    (
                        (allEdgeData[edgeI][0] > region0)
                     || (allEdgeData[edgeI][1] > region1)
                    )
                    {
                        allEdgeData[edgeI] = labelPair(region0, region1);
                        if (!isChangedEdge[edgeI])
                        {
                            changedEdges.append(edgeI);
                            isChangedEdge[edgeI] = true;
                        }
                    }
                }
            }
        }

        syncEdges
        (
            globalData,

            patchEdges,
            coupledEdges,
            sameEdgeOrientation,
            syncNonCollocated,

            isChangedEdge,
            changedEdges,
            allEdgeData
        );


        label nChangedEdges = returnReduce(changedEdges.size(), sumOp<label>());
        if (nChangedEdges == 0)
        {
            break;
        }
    }



    // Assign local regions
    // ~~~~~~~~~~~~~~~~~~~~

    // Calculate addressing from global region back to local region
    pointLocalRegions.setSize(patch.size());
    Map<label> globalToLocalRegion(globalRegions.localSize()/4);
    DynamicList<label> dynLocalToGlobalRegion(globalToLocalRegion.size());
    forAll(patch.localFaces(), facei)
    {
        const face& f = patch.localFaces()[facei];
        face& pRegions = pointLocalRegions[facei];
        pRegions.setSize(f.size());
        forAll(f, fp)
        {
            label globalRegionI = pointGlobalRegions[facei][fp];

            Map<label>::iterator fnd = globalToLocalRegion.find(globalRegionI);

            if (fnd != globalToLocalRegion.end())
            {
                // Already encountered this global region. Assign same local one
                pRegions[fp] = fnd();
            }
            else
            {
                // Region not yet seen. Create new one
                label localRegionI = globalToLocalRegion.size();
                pRegions[fp] = localRegionI;
                globalToLocalRegion.insert(globalRegionI, localRegionI);
                dynLocalToGlobalRegion.append(globalRegionI);
            }
        }
    }
    localToGlobalRegion.transfer(dynLocalToGlobalRegion);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::createShellMesh::createShellMesh
(
    const primitiveFacePatch& patch,
    const faceList& pointRegions,
    const labelList& regionPoints
)
:
    patch_(patch),
    pointRegions_(pointRegions),
    regionPoints_(regionPoints)
{
    if (pointRegions_.size() != patch_.size())
    {
        FatalErrorInFunction
            << "nFaces:" << patch_.size()
            << " pointRegions:" << pointRegions.size()
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::createShellMesh::setRefinement
(
    const pointField& firstLayerDisp,
    const scalar expansionRatio,
    const label nLayers,
    const labelList& topPatchID,
    const labelList& bottomPatchID,
    const labelListList& extrudeEdgePatches,
    polyTopoChange& meshMod
)
{
    if (firstLayerDisp.size() != regionPoints_.size())
    {
        FatalErrorInFunction
            << "nRegions:" << regionPoints_.size()
            << " firstLayerDisp:" << firstLayerDisp.size()
            << exit(FatalError);
    }

    if
    (
        topPatchID.size() != patch_.size()
     && bottomPatchID.size() != patch_.size()
    )
    {
        FatalErrorInFunction
            << "nFaces:" << patch_.size()
            << " topPatchID:" << topPatchID.size()
            << " bottomPatchID:" << bottomPatchID.size()
            << exit(FatalError);
    }

    if (extrudeEdgePatches.size() != patch_.nEdges())
    {
        FatalErrorInFunction
            << "nEdges:" << patch_.nEdges()
            << " extrudeEdgePatches:" << extrudeEdgePatches.size()
            << exit(FatalError);
    }



    // From cell to patch (trivial)
    DynamicList<label> cellToFaceMap(nLayers*patch_.size());
    // From face to patch+turning index
    DynamicList<label> faceToFaceMap
    (
        (nLayers+1)*(patch_.size()+patch_.nEdges())
    );
    // From face to patch edge index
    DynamicList<label> faceToEdgeMap(nLayers*(patch_.nEdges()+patch_.nEdges()));
    // From point to patch point index
    DynamicList<label> pointToPointMap((nLayers+1)*patch_.nPoints());


    // Introduce new cell for every face
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList addedCells(nLayers*patch_.size());
    forAll(patch_, facei)
    {
        for (label layerI = 0; layerI < nLayers; layerI++)
        {
            addedCells[nLayers*facei+layerI] = meshMod.addCell
            (
                -1,                     // masterPointID
                -1,                     // masterEdgeID
                -1,                     // masterFaceID
                cellToFaceMap.size(),   // masterCellID
                -1                      // zoneID
            );
            cellToFaceMap.append(facei);
        }
    }


    // Introduce original points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    // Original point numbers in local point ordering so no need to store.
    forAll(patch_.localPoints(), pointi)
    {
        // label addedPointi =
        meshMod.addPoint
        (
            patch_.localPoints()[pointi],   // point
            pointToPointMap.size(),         // masterPointID
            -1,                             // zoneID
            true                            // inCell
        );
        pointToPointMap.append(pointi);

        // Pout<< "Added bottom point " << addedPointi
        //    << " at " << patch_.localPoints()[pointi]
        //    << "  from point " << pointi
        //    << endl;
    }


    // Introduce new points (one for every region)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList addedPoints(nLayers*regionPoints_.size());
    forAll(regionPoints_, regionI)
    {
        label pointi = regionPoints_[regionI];

        point pt = patch_.localPoints()[pointi];
        point disp = firstLayerDisp[regionI];
        for (label layerI = 0; layerI < nLayers; layerI++)
        {
            pt += disp;

            addedPoints[nLayers*regionI+layerI] = meshMod.addPoint
            (
                pt,                     // point
                pointToPointMap.size(), // masterPointID - used only addressing
                -1,                     // zoneID
                true                    // inCell
            );
            pointToPointMap.append(pointi);

            disp *= expansionRatio;
        }
    }


    // Add face on bottom side
    forAll(patch_.localFaces(), facei)
    {
        meshMod.addFace
        (
            patch_.localFaces()[facei].reverseFace(),// vertices
            addedCells[nLayers*facei],  // own
            -1,                         // nei
            -1,                         // masterPointID
            -1,                         // masterEdgeID
            faceToFaceMap.size(),       // masterFaceID : current facei
            true,                       // flipFaceFlux
            bottomPatchID[facei],       // patchID
            -1,                         // zoneID
            false                       // zoneFlip
        );
        faceToFaceMap.append(-facei-1); // points to flipped original face
        faceToEdgeMap.append(-1);

        // const face newF(patch_.localFaces()[facei].reverseFace());
        // Pout<< "Added bottom face "
        //    << newF
        //    << " coords:" << UIndirectList<point>(meshMod.points(), newF)
        //    << " own " << addedCells[facei]
        //    << " patch:" << bottomPatchID[facei]
        //    << "  at " << patch_.faceCentres()[facei]
        //    << endl;
    }

    // Add in between faces and face on top
    forAll(patch_.localFaces(), facei)
    {
        // Get face in original ordering
        const face& f = patch_.localFaces()[facei];

        face newF(f.size());

        for (label layerI = 0; layerI < nLayers; layerI++)
        {
            // Pick up point based on region and layer
            forAll(f, fp)
            {
                label region = pointRegions_[facei][fp];
                newF[fp] = addedPoints[region*nLayers+layerI];
            }

            label own = addedCells[facei*nLayers+layerI];
            label nei;
            label patchi;
            if (layerI == nLayers-1)
            {
                nei = -1;
                patchi = topPatchID[facei];
            }
            else
            {
                nei = addedCells[facei*nLayers+layerI+1];
                patchi = -1;
            }

            meshMod.addFace
            (
                newF,                       // vertices
                own,                        // own
                nei,                        // nei
                -1,                         // masterPointID
                -1,                         // masterEdgeID
                faceToFaceMap.size(),       // masterFaceID : current facei
                false,                      // flipFaceFlux
                patchi,                     // patchID
                -1,                         // zoneID
                false                       // zoneFlip
            );
            faceToFaceMap.append(facei+1);  // unflipped
            faceToEdgeMap.append(-1);

            // Pout<< "Added in between face " << newF
            //    << " coords:" << UIndirectList<point>(meshMod.points(), newF)
            //    << " at layer " << layerI
            //    << " own " << own
            //    << " nei " << nei
            //    << "  at " << patch_.faceCentres()[facei]
            //    << endl;
        }
    }


    // Add side faces
    // ~~~~~~~~~~~~~~
    // Note that we loop over edges multiple times so for edges with
    // two cyclic faces they get added in two passes (for correct ordering)

    // Pass1. Internal edges and first face of other edges
    forAll(extrudeEdgePatches, edgeI)
    {
        const labelList& eFaces = patch_.edgeFaces()[edgeI];
        const labelList& ePatches = extrudeEdgePatches[edgeI];

        if (ePatches.size() == 0)
        {
            // Internal face
            if (eFaces.size() != 2)
            {
                FatalErrorInFunction
                    << "edge:" << edgeI
                    << " not internal but does not have side-patches defined."
                    << exit(FatalError);
            }
        }
        else
        {
            if (eFaces.size() != ePatches.size())
            {
                FatalErrorInFunction
                    << "external/feature edge:" << edgeI
                    << " has " << eFaces.size() << " connected extruded faces "
                    << " but only " << ePatches.size()
                    << " boundary faces defined." << exit(FatalError);
            }
        }



        // Make face pointing in to eFaces[0] so out of new master face
        const face& f = patch_.localFaces()[eFaces[0]];
        const edge& e = patch_.edges()[edgeI];

        label fp0 = findIndex(f, e[0]);
        label fp1 = f.fcIndex(fp0);

        if (f[fp1] != e[1])
        {
            fp1 = fp0;
            fp0 = f.rcIndex(fp1);
        }

        face newF(4);

        for (label layerI = 0; layerI < nLayers; layerI++)
        {
            label region0 = pointRegions_[eFaces[0]][fp0];
            label region1 = pointRegions_[eFaces[0]][fp1];

            // Pick up points with correct normal
            if (layerI == 0)
            {
                newF[0] = f[fp0];
                newF[1] = f[fp1];
                newF[2] = addedPoints[nLayers*region1+layerI];
                newF[3] = addedPoints[nLayers*region0+layerI];
            }
            else
            {
                newF[0] = addedPoints[nLayers*region0+layerI-1];
                newF[1] = addedPoints[nLayers*region1+layerI-1];
                newF[2] = addedPoints[nLayers*region1+layerI];
                newF[3] = addedPoints[nLayers*region0+layerI];
            }

            // Optionally rotate so e[0] is always 0th vertex. Note that
            // this normally is automatically done by coupled face ordering
            // but with NOORDERING we have to do it ourselves.
            if (f[fp0] != e[0])
            {
                // rotate one back to get newF[1] (originating from e[0])
                // into newF[0]
                label v0 = newF[0];
                for (label i = 0; i < newF.size()-1; i++)
                {
                    newF[i] = newF[newF.fcIndex(i)];
                }
                newF.last() = v0;
            }


            label minCelli = addedCells[nLayers*eFaces[0]+layerI];
            label maxCelli;
            label patchi;
            if (ePatches.size() == 0)
            {
                maxCelli = addedCells[nLayers*eFaces[1]+layerI];
                if (minCelli > maxCelli)
                {
                    // Swap
                    Swap(minCelli, maxCelli);
                    newF = newF.reverseFace();
                }
                patchi = -1;
            }
            else
            {
                maxCelli = -1;
                patchi = ePatches[0];
            }

            //{
            //    Pout<< "Adding from face:" << patch_.faceCentres()[eFaces[0]]
            //        << " from edge:"
            //        << patch_.localPoints()[f[fp0]]
            //        << patch_.localPoints()[f[fp1]]
            //        << " at layer:" << layerI
            //        << " with new points:" << newF
            //        << " locations:"
            //        << UIndirectList<point>(meshMod.points(), newF)
            //        << " own:" << minCelli
            //        << " nei:" << maxCelli
            //        << endl;
            //}


            // newF already outwards pointing.
            meshMod.addFace
            (
                newF,                   // vertices
                minCelli,               // own
                maxCelli,               // nei
                -1,                     // masterPointID
                -1,                     // masterEdgeID
                faceToFaceMap.size(),   // masterFaceID
                false,                  // flipFaceFlux
                patchi,                 // patchID
                -1,                     // zoneID
                false                   // zoneFlip
            );
            faceToFaceMap.append(0);
            faceToEdgeMap.append(edgeI);
        }
    }

    // Pass2. Other faces of boundary edges
    forAll(extrudeEdgePatches, edgeI)
    {
        const labelList& eFaces = patch_.edgeFaces()[edgeI];
        const labelList& ePatches = extrudeEdgePatches[edgeI];

        if (ePatches.size() >= 2)
        {
            for (label i = 1; i < ePatches.size(); i++)
            {
                // Extrude eFaces[i]
                label minFacei = eFaces[i];

                // Make face pointing in to eFaces[0] so out of new master face
                const face& f = patch_.localFaces()[minFacei];

                const edge& e = patch_.edges()[edgeI];
                label fp0 = findIndex(f, e[0]);
                label fp1 = f.fcIndex(fp0);

                if (f[fp1] != e[1])
                {
                    fp1 = fp0;
                    fp0 = f.rcIndex(fp1);
                }

                face newF(4);
                for (label layerI = 0; layerI < nLayers; layerI++)
                {
                    label region0 = pointRegions_[minFacei][fp0];
                    label region1 = pointRegions_[minFacei][fp1];

                    if (layerI == 0)
                    {
                        newF[0] = f[fp0];
                        newF[1] = f[fp1];
                        newF[2] = addedPoints[nLayers*region1+layerI];
                        newF[3] = addedPoints[nLayers*region0+layerI];
                    }
                    else
                    {
                        newF[0] = addedPoints[nLayers*region0+layerI-1];
                        newF[1] = addedPoints[nLayers*region1+layerI-1];
                        newF[2] = addedPoints[nLayers*region1+layerI];
                        newF[3] = addedPoints[nLayers*region0+layerI];
                    }


                    // Optionally rotate so e[0] is always 0th vertex. Note that
                    // this normally is automatically done by coupled face
                    // ordering but with NOORDERING we have to do it ourselves.
                    if (f[fp0] != e[0])
                    {
                        // rotate one back to get newF[1] (originating
                        // from e[0]) into newF[0].
                        label v0 = newF[0];
                        for (label i = 0; i < newF.size()-1; i++)
                        {
                            newF[i] = newF[newF.fcIndex(i)];
                        }
                        newF.last() = v0;
                    }
                    ////if (ePatches.size() == 0)
                    //{
                    //    Pout<< "Adding from MULTI face:"
                    //        << patch_.faceCentres()[minFacei]
                    //        << " from edge:"
                    //        << patch_.localPoints()[f[fp0]]
                    //        << patch_.localPoints()[f[fp1]]
                    //        << " at layer:" << layerI
                    //        << " with new points:" << newF
                    //        << " locations:"
                    //        << UIndirectList<point>(meshMod.points(), newF)
                    //        << endl;
                    //}

                    // newF already outwards pointing.
                    meshMod.addFace
                    (
                        newF,                   // vertices
                        addedCells[nLayers*minFacei+layerI],   // own
                        -1,                     // nei
                        -1,                     // masterPointID
                        -1,                     // masterEdgeID
                        faceToFaceMap.size(),   // masterFaceID
                        false,                  // flipFaceFlux
                        ePatches[i],            // patchID
                        -1,                     // zoneID
                        false                   // zoneFlip
                    );
                    faceToFaceMap.append(0);
                    faceToEdgeMap.append(edgeI);
                }
            }
        }
    }


    cellToFaceMap_.transfer(cellToFaceMap);
    faceToFaceMap_.transfer(faceToFaceMap);
    faceToEdgeMap_.transfer(faceToEdgeMap);
    pointToPointMap_.transfer(pointToPointMap);
}


void Foam::createShellMesh::updateMesh(const mapPolyMesh& map)
{
    inplaceReorder(map.reverseCellMap(), cellToFaceMap_);
    inplaceReorder(map.reverseFaceMap(), faceToFaceMap_);
    inplaceReorder(map.reverseFaceMap(), faceToEdgeMap_);
    inplaceReorder(map.reversePointMap(), pointToPointMap_);
}


// ************************************************************************* //
