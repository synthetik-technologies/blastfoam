/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify i
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

#include "cellVolumeWeightCellCellStencil.H"
#include "addToRunTimeSelectionTable.H"
#include "OBJstream.H"
#include "Time.H"
#include "oversetMeshToMesh.H"
#include "oversetCellVolumeWeightMethod.H"
#include "fvMeshSubset.H"
#include "regionSplit.H"
#include "globalIndex.H"
#include "oversetFvPatch.H"
#include "zeroGradientFvPatchFields.H"
#include "syncTools.H"
#include "dynamicOversetBlastFvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cellCellStencils
{
    defineTypeNameAndDebug(cellVolumeWeight, 0);
    addToRunTimeSelectionTable(cellCellStencil, cellVolumeWeight, mesh);
}
}

Foam::scalar
Foam::cellCellStencils::cellVolumeWeight::defaultOverlapTolerance_ = 1e-6;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cellCellStencils::cellVolumeWeight::walkFront
(
    const scalar layerRelax,
    labelList& allCellTypes,
    scalarField& allWeight
) const
{
    // Current front
    PackedBoolList isFront(mesh_.nFaces());
    // unused: bitSet doneCell(mesh_.nCells());

    const fvBoundaryMesh& fvm = mesh_.boundary();


    // 'overset' patches

    forAll(fvm, patchI)
    {
        if (isA<oversetFvPatch>(fvm[patchI]))
        {
            //Pout<< "Storing faces on patch " << fvm[patchI].name() << endl;
            forAll(fvm[patchI], i)
            {
                isFront.set(fvm[patchI].start()+i);
            }
        }
    }


    // Outside of 'hole' region
    {
        const labelList& own = mesh_.faceOwner();
        const labelList& nei = mesh_.faceNeighbour();

        for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
        {
            label ownType = allCellTypes[own[faceI]];
            label neiType = allCellTypes[nei[faceI]];
            if
            (
                 (ownType == HOLE && neiType != HOLE)
              || (ownType != HOLE && neiType == HOLE)
            )
            {
                //Pout<< "Front at face:" << faceI
                //    << " at:" << mesh_.faceCentres()[faceI] << endl;
                isFront.set(faceI);
            }
        }

        labelList nbrCellTypes;
        syncTools::swapBoundaryCellList(mesh_, allCellTypes, nbrCellTypes);

        for
        (
            label faceI = mesh_.nInternalFaces();
            faceI < mesh_.nFaces();
            faceI++
        )
        {
            label ownType = allCellTypes[own[faceI]];
            label neiType = nbrCellTypes[faceI-mesh_.nInternalFaces()];

            if
            (
                 (ownType == HOLE && neiType != HOLE)
              || (ownType != HOLE && neiType == HOLE)
            )
            {
                //Pout<< "Front at coupled face:" << faceI
                //    << " at:" << mesh_.faceCentres()[faceI] << endl;
                isFront.set(faceI);
            }
        }
    }



    // Current interpolation fraction
    scalar fraction = 1.0;

    while (fraction > SMALL && returnReduce(isFront.count(), sumOp<label>()))
    {
        // Interpolate cells on front

        Info<< "Front : fraction:" << fraction
            << " size:" << returnReduce(isFront.count(), sumOp<label>())
            << endl;

        PackedBoolList newIsFront(mesh_.nFaces());
        forAll(isFront, faceI)
        {
            if (isFront.get(faceI))
            {
                label own = mesh_.faceOwner()[faceI];
                if (allCellTypes[own] != HOLE)
                {
                    if (allWeight[own] < fraction)
                    {
                        allWeight[own] = fraction;

                        //if (debug)
                        //{
                        //    Pout<< "    setting cell "
                        //        << mesh_.cellCentres()[own]
                        //        << " to " << fraction << endl;
                        //}
                        allCellTypes[own] = INTERPOLATED;
                        newIsFront.set(mesh_.cells()[own]);
                    }
                }
                if (mesh_.isInternalFace(faceI))
                {
                    label nei = mesh_.faceNeighbour()[faceI];
                    if (allCellTypes[nei] != HOLE)
                    {
                        if (allWeight[nei] < fraction)
                        {
                            allWeight[nei] = fraction;

                            //if (debug)
                            //{
                            //    Pout<< "    setting cell "
                            //        << mesh_.cellCentres()[nei]
                            //        << " to " << fraction << endl;
                            //}

                            allCellTypes[nei] = INTERPOLATED;
                            newIsFront.set(mesh_.cells()[nei]);
                        }
                    }
                }
            }
        }

        syncTools::syncFaceList(mesh_, newIsFront, orEqOp<unsigned int>());

        isFront.transfer(newIsFront);

        fraction -= layerRelax;
    }
}


void Foam::cellCellStencils::cellVolumeWeight::findHoles
(
    const globalIndex& globalCells,
    const fvMesh& mesh,
    const labelList& zoneID,
    const labelListList& stencil,
    labelList& cellTypes
) const
{
    const fvBoundaryMesh& pbm = mesh.boundary();
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();


    // The input cellTypes will be
    // - HOLE           : cell part covered by other-mesh patch
    // - INTERPOLATED   : cell fully covered by other-mesh patch
    //                    or next to 'overset' patch
    // - CALCULATED     : otherwise
    //
    // so we start a walk from our patches and any cell we cannot reach
    // (because we walk is stopped by other-mesh patch) is a hole.


    boolList isBlockedFace(mesh.nFaces(), false);
    for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        label ownType = cellTypes[own[faceI]];
        label neiType = cellTypes[nei[faceI]];
        if
        (
             (ownType == HOLE && neiType != HOLE)
          || (ownType != HOLE && neiType == HOLE)
        )
        {
            isBlockedFace[faceI] = true;
        }
    }

    labelList nbrCellTypes;
    syncTools::swapBoundaryCellList(mesh, cellTypes, nbrCellTypes);

    for (label faceI = mesh.nInternalFaces(); faceI < mesh.nFaces(); faceI++)
    {
        label ownType = cellTypes[own[faceI]];
        label neiType = nbrCellTypes[faceI-mesh.nInternalFaces()];

        if
        (
             (ownType == HOLE && neiType != HOLE)
          || (ownType != HOLE && neiType == HOLE)
        )
        {
            isBlockedFace[faceI] = true;
        }
    }

    regionSplit cellRegion(mesh, isBlockedFace);

    Info<< typeName << " : detected " << cellRegion.nRegions()
        << " mesh regions after overset" << nl << endl;



    // Now we'll have a mesh split according to where there are cells
    // covered by the other-side patches. See what we can reach from our
    // real patches

    //  0 : region not yet determined
    //  1 : borders blockage so is not ok (but can be overridden by real
    //      patch)
    //  2 : has real patch in it so is reachable
    labelList regionType(cellRegion.nRegions(), Zero);


    // See if any regions borders blockage. Note: isBlockedFace is already
    // parallel synchronised.
    {
        for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
        {
            if (isBlockedFace[faceI])
            {
                label ownRegion = cellRegion[own[faceI]];

                if (cellTypes[own[faceI]] != HOLE)
                {
                    if (regionType[ownRegion] == 0)
                    {
                        //Pout<< "Mark region:" << ownRegion
                        //    << " on zone:" << zoneID[own[faceI]]
                        //    << " as next to blockage at:"
                        //    << mesh.faceCentres()[faceI] << endl;
                        regionType[ownRegion] = 1;
                    }
                }

                label neiRegion = cellRegion[nei[faceI]];

                if (cellTypes[nei[faceI]] != HOLE)
                {
                    if (regionType[neiRegion] == 0)
                    {
                        //Pout<< "Mark region:" << neiRegion
                        //    << " on zone:" << zoneID[nei[faceI]]
                        //    << " as next to blockage at:"
                        //    << mesh.faceCentres()[faceI] << endl;
                        regionType[neiRegion] = 1;
                    }
                }
            }
        }
        for
        (
            label faceI = mesh.nInternalFaces();
            faceI < mesh.nFaces();
            faceI++
        )
        {
            if (isBlockedFace[faceI])
            {
                label ownRegion = cellRegion[own[faceI]];

                if (regionType[ownRegion] == 0)
                {
                    //Pout<< "Mark region:" << ownRegion
                    //    << " on zone:" << zoneID[own[faceI]]
                    //    << " as next to blockage at:"
                    //    << mesh.faceCentres()[faceI] << endl;
                    regionType[ownRegion] = 1;
                }
            }
        }
    }


    // Override with real patches
    forAll(pbm, patchI)
    {
        const fvPatch& fvp = pbm[patchI];

        if (isA<oversetFvPatch>(fvp))
        {}
        else if (!fvPatch::constraintType(fvp.type()))
        {
            //Pout<< "Proper patch " << fvp.name() << " of type " << fvp.type()
            //    << endl;

            const labelList& fc = fvp.faceCells();
            forAll(fc, i)
            {
                label regionI = cellRegion[fc[i]];

                if (cellTypes[fc[i]] != HOLE && regionType[regionI] != 2)
                {
                    //Pout<< "reachable region : " << regionI
                    //    << " at cell " << mesh.cellCentres()[fc[i]]
                    //    << " on zone " << zoneID[fc[i]] << endl;
                    regionType[regionI] = 2;
                }
            }
        }
    }

    // Now we've handled
    // - cells next to blocked cells
    // - coupled boundaries
    // Only thing to handle is the interpolation between regions


    labelListList compactStencil(stencil);
    List<Map<label>> compactMap;
    mapDistribute map(globalCells, compactStencil, compactMap);

    while (true)
    {
        // Synchronise region status on processors
        // (could instead swap status through processor patches)
        Pstream::listCombineGather(regionType, maxEqOp<label>());
        Pstream::listCombineScatter(regionType);

        // Communicate region status through interpolative cells
        labelList cellRegionType(UIndirectList<label>(regionType, cellRegion));
        map.distribute(cellRegionType);


        label nChanged = 0;
        forAll(pbm, patchI)
        {
            const fvPatch& fvp = pbm[patchI];

            if (isA<oversetFvPatch>(fvp))
            {
                const labelUList& fc = fvp.faceCells();
                forAll(fc, i)
                {
                    label cellI = fc[i];
                    label regionI = cellRegion[cellI];

                    if (regionType[regionI] != 2)
                    {
                        const labelList& slots = compactStencil[cellI];
                        forAll(slots, i)
                        {
                            label otherType = cellRegionType[slots[i]];

                            if (otherType == 2)
                            {
                                //Pout<< "Reachable through interpolation : "
                                //    << regionI << " at cell "
                                //    << mesh.cellCentres()[cellI] << endl;
                                regionType[regionI] = 2;
                                nChanged++;
                                break;
                            }
                        }
                    }
                }
            }
        }


        reduce(nChanged, sumOp<label>());
        if (nChanged == 0)
        {
            break;
        }
    }


    // See which regions have not been visited (regionType == 1)
    forAll(cellRegion, cellI)
    {
        label type = regionType[cellRegion[cellI]];
        if (type == 1 && cellTypes[cellI] != HOLE)
        {
            cellTypes[cellI] = HOLE;
        }
    }
}


void Foam::cellCellStencils::cellVolumeWeight::markPatchCells
(
    const fvMesh& mesh,
    const labelList& cellMap,
    labelList& patchCellTypes
) const
{
    const fvBoundaryMesh& pbm = mesh.boundary();

    forAll(pbm, patchI)
    {
        const fvPatch& fvp = pbm[patchI];
        const labelList& fc = fvp.faceCells();

        if (isA<oversetFvPatch>(fvp))
        {
            //Pout<< "Marking cells on overset patch " << fvp.name() << endl;
            forAll(fc, i)
            {
                label cellI = fc[i];
                patchCellTypes[cellMap[cellI]] = OVERSET;
            }
        }
        else if (!fvPatch::constraintType(fvp.type()))
        {
            //Pout<< "Marking cells on proper patch " << fvp.name()
            //    << " with type " << fvp.type() << endl;
            forAll(fc, i)
            {
                label cellI = fc[i];
                if (patchCellTypes[cellMap[cellI]] != OVERSET)
                {
                    patchCellTypes[cellMap[cellI]] = PATCH;
                }
            }
        }
    }
}


void Foam::cellCellStencils::cellVolumeWeight::interpolatePatchTypes
(
    const labelListList& addressing,
    const labelList& patchTypes,
    labelList& result
) const
{
    forAll(result, cellI)
    {
        const labelList& slots = addressing[cellI];
        forAll(slots, i)
        {
            label type = patchTypes[slots[i]];

            if (type == OVERSET)
            {
                // 'overset' overrides anything
                result[cellI] = OVERSET;
                break;
            }
            else if (type == PATCH)
            {
                // 'patch' overrides -1 and 'other'
                result[cellI] = PATCH;
            }
            else if (result[cellI] == -1)
            {
                // 'other' overrides -1 only
                result[cellI] = OTHER;
            }
        }
    }
}


void Foam::cellCellStencils::cellVolumeWeight::interpolatePatchTypes
(
    const autoPtr<mapDistribute>& mapPtr,
    const labelListList& addressing,
    const labelList& patchTypes,
    labelList& result
) const
{
    if (result.size() != addressing.size())
    {
        FatalErrorInFunction << "result:" << result.size()
            << " addressing:" << addressing.size() << exit(FatalError);
    }


    // Initialise to not-mapped
    result = -1;

    if (mapPtr.valid())
    {
        // Pull remote data into order of addressing
        labelList work(patchTypes);
        mapPtr().distribute(work);

        interpolatePatchTypes(addressing, work, result);
    }
    else
    {
        interpolatePatchTypes(addressing, patchTypes, result);
    }
}


void Foam::cellCellStencils::cellVolumeWeight::combineCellTypes
(
    const label subZoneID,
    const fvMesh& subMesh,
    const labelList& subCellMap,

    const label donorZoneID,
    const labelListList& addressing,
    const List<scalarList>& weights,
    const labelList& otherCells,
    const labelList& interpolatedOtherPatchTypes,

    labelListList& allStencil,
    scalarListList& allWeights,
    labelList& allCellTypes,
    labelList& allDonorID
) const
{
    forAll(subCellMap, subCellI)
    {
        label cellI = subCellMap[subCellI];

        bool validDonors = true;
        switch (interpolatedOtherPatchTypes[subCellI])
        {
            case -1:
            {
                validDonors = false;
            }
            break;

            case OTHER:
            {
                // No patch interaction so keep valid
            }
            break;

            case PATCH:
            {
                // Patch-patch interaction... For now disable always
                allCellTypes[cellI] = HOLE;
                validDonors = false;

                // Alternative is to look at the amount of overlap but this
                // is not very robust
                //if (allCellTypes[cellI] != HOLE)
                //{
                //    scalar overlapVol = sum(weights[subCellI]);
                //    scalar v = mesh_.V()[cellI];
                //    if (overlapVol < (1.0-overlapTolerance_)*v)
                //    {
                //        //Pout<< "** Patch overlap:" << cellI
                //        //    << " at:" << mesh_.cellCentres()[cellI] << endl;
                //        allCellTypes[cellI] = HOLE;
                //        validDonors = false;
                //    }
                //}
            }
            break;

            case OVERSET:
            {
                validDonors = false;
            }
            break;
        }


        if (validDonors)
        {
            // There are a few possible choices how to choose between multiple
            // donor candidates:
            // 1 highest overlap volume. However this is generally already
            //   99.9% so you're just measuring truncation error.
            // 2 smallest donors cells or most donor cells. This is quite
            //   often done but can cause switching of donor zone from one
            //   time step to the other if the donor meshes are non-uniform
            //   and the acceptor cells just happens to be sweeping through
            //   some small donor cells.
            // 3 nearest zoneID. So zone 0 preferentially interpolates from
            //   zone 1, zone 1 preferentially from zone 2 etc.

            //- Option 1:
            //scalar currentVol = sum(allWeights[cellI]);
            //if (overlapVol[subCellI] > currentVol)

            //- Option 3:
            label currentDiff = mag(subZoneID-allDonorID[cellI]);
            label thisDiff = mag(subZoneID-donorZoneID);

            if
            (
                allDonorID[cellI] == -1
             || (thisDiff < currentDiff)
             || (thisDiff == currentDiff && donorZoneID > allDonorID[cellI])
            )
            {
                allWeights[cellI] = weights[subCellI];
                allStencil[cellI] =
                    UIndirectList<label>(otherCells, addressing[subCellI]);
                allDonorID[cellI] = donorZoneID;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellCellStencils::cellVolumeWeight::cellVolumeWeight
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    cellCellStencil(mesh),
    dict_(dict),
    overlapTolerance_(defaultOverlapTolerance_),
    cellTypes_(labelList(mesh.nCells(), CALCULATED)),
    interpolationCells_(0),
    cellInterpolationMap_(),
    cellStencil_(0),
    cellInterpolationWeights_(0),
    cellInterpolationWeight_
    (
        IOobject
        (
            "cellInterpolationWeight",
            mesh_.facesInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    )
{
    // Protect local fields from interpolation
    nonInterpolatedFields_.insert("cellTypes");
    nonInterpolatedFields_.insert("cellInterpolationWeight");

    // For convenience also suppress frequently used displacement field
    nonInterpolatedFields_.insert("cellDisplacement");
    nonInterpolatedFields_.insert("grad(cellDisplacement)");
    const word w("snGradCorr(cellDisplacement)");
    const word d("((viscosity*faceDiffusivity)*magSf)");
    nonInterpolatedFields_.insert("surfaceIntegrate(("+d+"*"+w+"))");

    // Read zoneID
    this->zoneID();

    // Read old-time cellTypes
    IOobject io
    (
        "cellTypes",
        mesh_.time().timeName(),
        mesh_,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );
    if (io.typeHeaderOk<volScalarField>(true))
    {
        if (debug)
        {
            Pout<< "Reading cellTypes from time " << mesh_.time().timeName()
                << endl;
        }

        const volScalarField volCellTypes(io, mesh_);
        forAll(volCellTypes, celli)
        {
            // Round to integer
            cellTypes_[celli] = volCellTypes[celli];
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellCellStencils::cellVolumeWeight::~cellVolumeWeight()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cellCellStencils::cellVolumeWeight::updateMesh
(
    const mapPolyMesh& map
)
{
    // Map data
    const labelList& cellMap = map.cellMap();

    labelList newCellType(cellMap.size());
    forAll(cellMap, newCelli)
    {
        label oldCelli = cellMap[newCelli];

        if (oldCelli == -1)
        {
            newCellType[newCelli] = -1;
        }
        else
        {
            newCellType[newCelli] = cellTypes_[oldCelli];
        }
    }
    cellTypes_.transfer(newCellType);
    cellCellStencil::updateMesh(map);
}


bool Foam::cellCellStencils::cellVolumeWeight::update() const
{
    scalar layerRelax(dict_.lookupOrDefault("layerRelax", 1.0));
    const labelIOList& zoneID = this->zoneID();

    label nZones = gMax(zoneID)+1;
    labelList nCellsPerZone(nZones, Zero);
    forAll(zoneID, cellI)
    {
        nCellsPerZone[zoneID[cellI]]++;
    }
    Pstream::listCombineGather(nCellsPerZone, plusEqOp<label>());
    Pstream::listCombineScatter(nCellsPerZone);


    Info<< typeName << " : detected " << nZones
        << " mesh regions" << nl << endl;


    PtrList<fvMeshSubset> meshParts(nZones);

    Info<< incrIndent;
    forAll(meshParts, zonei)
    {
        Info<< indent<< "zone:" << zonei << " nCells:"
            << nCellsPerZone[zonei] << nl;

        meshParts.set
        (
            zonei,
            new fvMeshSubset(mesh_)
        );
        meshParts[zonei].setLargeCellSubset(zoneID, zonei);
    }
    Info<< decrIndent;



    // Current best guess for cells. Includes best stencil. Weights should
    // add up to volume.
    labelList allCellTypes(mesh_.nCells(), CALCULATED);
    labelList allPatchTypes(mesh_.nCells(), OTHER);
    labelListList allStencil(mesh_.nCells());
    scalarListList allWeights(mesh_.nCells());
    // zoneID of donor
    labelList allDonorID(mesh_.nCells(), -1);


    // Marking patch cells
    forAll(meshParts, partI)
    {
        const fvMesh& partMesh = meshParts[partI].subMesh();
        const labelList& partCellMap = meshParts[partI].cellMap();

        // Mark cells with
        // - overset boundary
        // - other, proper boundary
        // - other cells
        Info<< "Marking patch-cells on zone " << partI << endl;
        markPatchCells(partMesh, partCellMap, allPatchTypes);
    }


    labelList nCells(count(3, allPatchTypes));
    Info<< nl
        << "After patch analysis : nCells : "
        << returnReduce(allPatchTypes.size(), sumOp<label>()) << nl
        << incrIndent
        << indent << "other  : " << nCells[OTHER] << nl
        << indent << "patch  : " << nCells[PATCH] << nl
        << indent << "overset: " << nCells[OVERSET] << nl
        << decrIndent << endl;

    globalIndex globalCells(mesh_.nCells());


    for (label srcI = 0; srcI < meshParts.size()-1; srcI++)
    {
        const fvMesh& srcMesh = meshParts[srcI].subMesh();
        const labelList& srcCellMap = meshParts[srcI].cellMap();

        for (label tgtI = srcI+1; tgtI < meshParts.size(); tgtI++)
        {
            const fvMesh& tgtMesh = meshParts[tgtI].subMesh();
            const labelList& tgtCellMap = meshParts[tgtI].cellMap();

            oversetMeshToMesh mapper
            (
                srcMesh,
                tgtMesh,
                oversetMeshToMesh::interpolationMethod::imCellVolumeWeight,
                HashTable<word>(0),     // patchMap,
                wordList(0),            // cuttingPatches
                oversetMeshToMesh::procMapMethod::pmAABB,
                false                   // do not normalise
            );


            {
                // Get tgt patch types on src mesh
                labelList interpolatedTgtPatchTypes(srcMesh.nCells(), -1);
                interpolatePatchTypes
                (
                    mapper.tgtMap(),            // How to get remote data local
                    mapper.srcToTgtCellAddr(),
                    labelList(UIndirectList<label>(allPatchTypes, tgtCellMap)),
                    interpolatedTgtPatchTypes
                );

                // Get target cell labels in global cell indexing (on overall
                // mesh)
                labelList tgtGlobalCells(tgtMesh.nCells());
                {
                    forAll(tgtCellMap, tgtCellI)
                    {
                        label cellI = tgtCellMap[tgtCellI];
                        tgtGlobalCells[tgtCellI] = globalCells.toGlobal(cellI);
                    }
                    if (mapper.tgtMap().valid())
                    {
                        mapper.tgtMap()->distribute(tgtGlobalCells);
                    }
                }
                combineCellTypes
                (
                    srcI,
                    srcMesh,
                    srcCellMap,

                    tgtI,
                    mapper.srcToTgtCellAddr(),
                    mapper.srcToTgtCellWght(),
                    tgtGlobalCells,
                    interpolatedTgtPatchTypes,

                    // Overall mesh data
                    allStencil,
                    allWeights,
                    allCellTypes,
                    allDonorID
                );
            }

            {
                // Get src patch types on tgt mesh
                labelList interpolatedSrcPatchTypes(tgtMesh.nCells(), -1);
                interpolatePatchTypes
                (
                    mapper.srcMap(),            // How to get remote data local
                    mapper.tgtToSrcCellAddr(),
                    labelList(UIndirectList<label>(allPatchTypes, srcCellMap)),
                    interpolatedSrcPatchTypes
                );

                labelList srcGlobalCells(srcMesh.nCells());
                {
                    forAll(srcCellMap, srcCellI)
                    {
                        label cellI = srcCellMap[srcCellI];
                        srcGlobalCells[srcCellI] = globalCells.toGlobal(cellI);
                    }
                    if (mapper.srcMap().valid())
                    {
                        mapper.srcMap()->distribute(srcGlobalCells);
                    }
                }

                combineCellTypes
                (
                    tgtI,
                    tgtMesh,
                    tgtCellMap,

                    srcI,
                    mapper.tgtToSrcCellAddr(),
                    mapper.tgtToSrcCellWght(),
                    srcGlobalCells,
                    interpolatedSrcPatchTypes,

                    // Overall mesh data
                    allStencil,
                    allWeights,
                    allCellTypes,
                    allDonorID
                );
            }
        }
    }


    if (debug)
    {
        tmp<volScalarField> tfld
        (
            createField(mesh_, "allCellTypes", allCellTypes)
        );
        tfld().write();
    }


    // Use the patch types and weights to decide what to do
    forAll(allPatchTypes, cellI)
    {
        if (allCellTypes[cellI] != HOLE)
        {
            switch (allPatchTypes[cellI])
            {
                case OVERSET:
                {
                    // Interpolate. Check if enough overlap
                    scalar v = mesh_.V()[cellI];
                    scalar overlapVol = sum(allWeights[cellI]);
                    if (overlapVol > overlapTolerance_*v)
                    {
                        allCellTypes[cellI] = INTERPOLATED;
                    }
                    else
                    {
                        //Pout<< "Holeing interpolated cell:" << cellI
                        //    << " at:" << mesh_.cellCentres()[cellI] << endl;
                        allCellTypes[cellI] = HOLE;
                        allWeights[cellI].clear();
                        allStencil[cellI].clear();
                    }
                    break;
                }
            }
        }
    }

/*
    // Knock out cell with insufficient interpolation weights
    forAll(allCellTypes, cellI)
    {
        if (allCellTypes[cellI] == INTERPOLATED)
        {
            scalar v = mesh_.V()[cellI];
            scalar overlapVol = sum(allWeights[cellI]);
            if (overlapVol < (1.0-overlapTolerance_)*v)
            {
                //Pout<< "Holeing cell:" << cellI
                //    << " at:" << mesh_.cellCentres()[cellI] << endl;
                allCellTypes[cellI] = HOLE;
                allWeights[cellI].clear();
                allStencil[cellI].clear();
            }
        }
    }
*/
    if (debug)
    {
        tmp<volScalarField> tfld
        (
            createField(mesh_, "allCellTypes_patch", allCellTypes)
        );
        //tfld.ref().correctBoundaryConditions();
        tfld().write();
    }


    // Check previous iteration cellTypes_ for any hole->calculated changes
    // If so set the cell either to interpolated (if there are donors) or
    // holes (if there are no donors). Note that any interpolated cell might
    // still be overwritten by the flood filling
    {
        label nCalculated = 0;

        forAll(cellTypes_, celli)
        {
            if (allCellTypes[celli] == CALCULATED && cellTypes_[celli] == HOLE)
            {
                if (allStencil[celli].size() == 0)
                {
                    // Reset to hole
                    allCellTypes[celli] = HOLE;
                    allWeights[celli].clear();
                    allStencil[celli].clear();
                }
                else
                {
                    allCellTypes[celli] = INTERPOLATED;
                    nCalculated++;
                }
            }
        }

        if (debug)
        {
            Pout<< "Detected " << nCalculated << " cells changing from hole"
                << " to calculated. Changed these to interpolated"
                << endl;
        }
    }


    // Mark unreachable bits
    findHoles(globalCells, mesh_, zoneID, allStencil, allCellTypes);

    if (debug)
    {
        tmp<volScalarField> tfld
        (
            createField(mesh_, "allCellTypes_hole", allCellTypes)
        );
        //tfld.ref().correctBoundaryConditions();
        tfld().write();
    }


    // Add buffer interpolation layer around holes
    scalarField allWeight(mesh_.nCells(), Zero);
    walkFront(layerRelax, allCellTypes, allWeight);

    if (debug)
    {
        tmp<volScalarField> tfld
        (
            createField(mesh_, "allCellTypes_front", allCellTypes)
        );
        //tfld.ref().correctBoundaryConditions();
        tfld().write();
    }

    // Normalise weights, Clear storage
    forAll(allCellTypes, cellI)
    {
        if (allCellTypes[cellI] == INTERPOLATED)
        {
            if (allWeight[cellI] < SMALL || allStencil[cellI].size() == 0)
            {
                //Pout<< "Clearing cell:" << cellI
                //    << " at:" << mesh_.cellCentres()[cellI] << endl;
                allWeights[cellI].clear();
                allStencil[cellI].clear();
                allWeight[cellI] = 0.0;
            }
            else
            {
                scalar s = sum(allWeights[cellI]);
                forAll(allWeights[cellI], i)
                {
                    allWeights[cellI][i] /= s;
                }
            }
        }
        else
        {
            allWeights[cellI].clear();
            allStencil[cellI].clear();
        }
    }


    // Write to volField for debugging
    if (debug)
    {
        volScalarField patchTypes
        (
            IOobject
            (
                "patchTypes",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(dimless, Zero),
            zeroGradientFvPatchScalarField::typeName
        );

        forAll(patchTypes.internalField(), cellI)
        {
            patchTypes[cellI] = allPatchTypes[cellI];
        }
        //patchTypes.correctBoundaryConditions();
        dynamicOversetBlastFvMesh::correctBoundaryConditions
        <
            volScalarField,
            oversetFvPatchField<scalar>
        >(patchTypes.boundaryFieldRef(), false);
        patchTypes.write();
    }
    if (debug)
    {
        volScalarField volTypes
        (
            IOobject
            (
                "cellTypes",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(dimless, Zero),
            zeroGradientFvPatchScalarField::typeName
        );

        forAll(volTypes.internalField(), cellI)
        {
            volTypes[cellI] = allCellTypes[cellI];
        }
        //volTypes.correctBoundaryConditions();
        dynamicOversetBlastFvMesh::correctBoundaryConditions
        <
            volScalarField,
            oversetFvPatchField<scalar>
        >(volTypes.boundaryFieldRef(), false);
        volTypes.write();
    }


//  // Check previous iteration cellTypes_ for any hole->calculated changes
//  {
//      label nCalculated = 0;
//
//      forAll(cellTypes_, celli)
//      {
//          if (allCellTypes[celli] == CALCULATED && cellTypes_[celli] == HOLE)
//          {
//              if (allStencil[celli].size() == 0)
//              {
//                  FatalErrorInFunction
//                      << "Cell:" << celli
//                      << " at:" << mesh_.cellCentres()[celli]
//                      << " zone:" << zoneID[celli]
//                      << " changed from hole to calculated"
//                      << " but there is no donor"
//                      << exit(FatalError);
//              }
//              else
//              {
//                  allCellTypes[celli] = INTERPOLATED;
//                  nCalculated++;
//              }
//          }
//      }
//
//      if (debug)
//      {
//          Pout<< "Detected " << nCalculated << " cells changing from hole"
//              << " to calculated. Changed these to interpolated"
//              << endl;
//      }
//  }


    cellTypes_.transfer(allCellTypes);
    cellStencil_.transfer(allStencil);
    cellInterpolationWeights_.transfer(allWeights);
    cellInterpolationWeight_.transfer(allWeight);
    //cellInterpolationWeight_.correctBoundaryConditions();
    dynamicOversetBlastFvMesh::correctBoundaryConditions
    <
        volScalarField,
        oversetFvPatchField<scalar>
    >(cellInterpolationWeight_.boundaryFieldRef(), false);

    DynamicList<label> interpolationCells;
    forAll(cellStencil_, cellI)
    {
        if (cellStencil_[cellI].size())
        {
            interpolationCells.append(cellI);
        }
    }
    interpolationCells_.transfer(interpolationCells);


    List<Map<label>> compactMap;
    cellInterpolationMap_.reset
    (
        new mapDistribute(globalCells, cellStencil_, compactMap)
    );

    // Dump interpolation stencil
    if (debug)
    {
        // Dump weight
        cellInterpolationWeight_.instance() = mesh_.time().timeName();
        cellInterpolationWeight_.write();


        mkDir(mesh_.time().timePath());
        OBJstream str(mesh_.time().timePath()/"stencil2.obj");
        Info<< typeName << " : dumping to " << str.name() << endl;
        pointField cc(mesh_.cellCentres());
        cellInterpolationMap().distribute(cc);

        forAll(interpolationCells_, compactI)
        {
            label cellI = interpolationCells_[compactI];
            const labelList& slots = cellStencil_[cellI];

            Pout<< "cellI:" << cellI << " at:"
                << mesh_.cellCentres()[cellI]
                << " calculated from slots:" << slots
                << " cc:" << UIndirectList<point>(cc, slots)
                << " weights:" << cellInterpolationWeights_[cellI]
                << endl;

            forAll(slots, i)
            {
                const point& donorCc = cc[slots[i]];
                const point& accCc = mesh_.cellCentres()[cellI];

                str.write(linePointRef(accCc, 0.1*accCc+0.9*donorCc));
            }
        }
    }


    {
        labelList nCells(count(3, cellTypes_));
        Info<< "Overset analysis : nCells : "
            << returnReduce(cellTypes_.size(), sumOp<label>()) << nl
            << incrIndent
            << indent << "calculated   : " << nCells[CALCULATED] << nl
            << indent << "interpolated : " << nCells[INTERPOLATED] << nl
            << indent << "hole         : " << nCells[HOLE] << nl
            << decrIndent << endl;
    }

    // Tbd: detect if anything changed. Most likely it did!
    return cellCellStencil::update();
}


void Foam::cellCellStencils::cellVolumeWeight::stencilWeights
(
    const point& sample,
    const pointList& donorCcs,
    scalarList& weights
) const
{
    // Inverse-distance weighting

    weights.setSize(donorCcs.size());
    scalar sum = 0.0;
    forAll(donorCcs, i)
    {
        scalar d = mag(sample-donorCcs[i]);

        if (d > ROOTVSMALL)
        {
            weights[i] = 1.0/d;
            sum += weights[i];
        }
        else
        {
            // Short circuit
            weights = 0.0;
            weights[i] = 1.0;
            return;
        }
    }
    forAll(weights, i)
    {
        weights[i] /= sum;
    }
}


// ************************************************************************* //
