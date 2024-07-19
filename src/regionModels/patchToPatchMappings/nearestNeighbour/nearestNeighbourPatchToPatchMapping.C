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

#include "nearestNeighbourPatchToPatchMapping.H"
#include "addToRunTimeSelectionTable.H"
#include "indexedOctree.H"
#include "treeDataPrimitivePatch.H"
#include "treeDataPoint.H"
#include "Time.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace patchToPatchMappings
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(nearestNeighbourPatchToPatchMapping, 0);
addToRunTimeSelectionTable
(
    patchToPatchMapping,
    nearestNeighbourPatchToPatchMapping,
    dictionary
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void nearestNeighbourPatchToPatchMapping::checkZoneSizes() const
{}


void nearestNeighbourPatchToPatchMapping::calcZoneAToZoneBFaceMap() const
{
    if (zoneAToZoneBFaceMapPtr_.valid())
    {
        FatalErrorInFunction
            << "List already set!" << abort(FatalError);
    }

    // Check zones are conformal
    checkZoneSizes();

    // Initialise map
    zoneAToZoneBFaceMapPtr_.set
    (
        new labelList(zoneB().size(), -1)
    );
    labelList& zoneToZoneMap = zoneAToZoneBFaceMapPtr_();

    // Perform N^2 search for corresponding faces
    // We will take 0.1% of the minEdgeLength as the exact match
    // relative tolerance
    treeBoundBox bbA(zoneA().localPoints());
    bbA = bbA.extend(1e-4);

    const scalar planarTol =
        indexedOctree<treeDataPrimitivePatch<standAlonePatch>>::
        perturbTol();
    indexedOctree<treeDataPrimitivePatch<standAlonePatch>> tree
    (
        treeDataPrimitivePatch<standAlonePatch>
        (
            false,
            zoneA(),
            planarTol
        ),
        bbA,
        10,
        10,
        3
    );

    scalar nds
    (
//         max
//         (
//             magSqr(boundBox(zoneB().localPoints()).span()),
//             magSqr(bbA.span())
//         )
        great
    );
    const vectorField& pCf = zoneB().faceCentres();
    forAll(pCf, facei)
    {
        const vector& pt = pCf[facei];
        pointIndexHit pIH = tree.findNearest(pt, nds);
        if (pIH.hit())
        {
            zoneToZoneMap[facei] = pIH.index();
        }
    }

    if (requireMatch_ && gMin(zoneToZoneMap) == -1)
    {
        FatalErrorInFunction
            << "Cannot calculate the map between interfaces!" << nl
            << "ZoneA and zoneB are not similar (zoneA patch ="
            << " " << globalPatchA().patchName() << ", zoneB "
            << "patch = " << globalPatchB().patchName() << ")" << endl
            << abort(FatalError);
    }
}


const labelList&
nearestNeighbourPatchToPatchMapping::zoneAToZoneBFaceMap() const
{
    if (zoneAToZoneBFaceMapPtr_.empty())
    {
        calcZoneAToZoneBFaceMap();
    }

    return zoneAToZoneBFaceMapPtr_;
}


void nearestNeighbourPatchToPatchMapping::calcZoneBToZoneAFaceMap() const
{
    if (zoneBToZoneAFaceMapPtr_.valid())
    {
        FatalErrorInFunction
            << "List already set!" << abort(FatalError);
    }

    // Check zones are conformal
    checkZoneSizes();

    // Initialise map
    zoneBToZoneAFaceMapPtr_.set
    (
        new labelList(labelList(zoneA().size(), -1))
    );
    labelList& zoneToZoneMap = zoneBToZoneAFaceMapPtr_();

    treeBoundBox bbB(zoneB().localPoints());
    bbB = bbB.extend(1e-4);

    const scalar planarTol =
        indexedOctree<treeDataPrimitivePatch<standAlonePatch>>::
        perturbTol();
    indexedOctree<treeDataPrimitivePatch<standAlonePatch>> tree
    (
        treeDataPrimitivePatch<standAlonePatch>
        (
            false,
            zoneB(),
            planarTol
        ),
        bbB,
        10,
        10,
        3
    );

    scalar nds
    (
//         max
//         (
//             magSqr(boundBox(zoneA().localPoints()).span()),
//             magSqr(bbB.span())
//         )
        great
    );
    const vectorField& pCf = zoneA().faceCentres();
    forAll(pCf, facei)
    {
        const vector& pt = pCf[facei];
        pointIndexHit pIH = tree.findNearest(pt, nds);
        if (pIH.hit())
        {
            zoneToZoneMap[facei] = pIH.index();
        }
    }

    if (requireMatch_ && gMin(zoneToZoneMap) == -1)
    {
        FatalErrorInFunction
            << "Cannot calculate the map between interfaces!" << nl
            << "ZoneA and zoneB are not similar (zoneA patch ="
            << " " << globalPatchA().patchName() << ", zoneB "
            << "patch = " << globalPatchB().patchName() << ")" << endl
            << abort(FatalError);
    }
}


const labelList&
nearestNeighbourPatchToPatchMapping::zoneBToZoneAFaceMap() const
{
    if (zoneBToZoneAFaceMapPtr_.empty())
    {
        calcZoneBToZoneAFaceMap();
    }

    return zoneBToZoneAFaceMapPtr_;
}


void nearestNeighbourPatchToPatchMapping::calcZoneAToZoneBPointMap() const
{
    if (zoneAToZoneBPointMapPtr_.valid())
    {
        FatalErrorInFunction
            << "List already set!" << abort(FatalError);
    }

    // Check zones are conformal
    checkZoneSizes();

    // Initialise map
    zoneAToZoneBPointMapPtr_.set
    (
        new labelList(labelList(zoneB().nPoints(), -1))
    );
    labelList& zoneToZoneMap = zoneAToZoneBPointMapPtr_();

    treeBoundBox bbA(zoneA().localPoints());
    bbA = bbA.extend(1e-4);

    indexedOctree<treeDataPoint> tree
    (
        treeDataPoint(zoneA().localPoints()),
        bbA,
        10,
        10,
        3
    );

    const scalar nds
    (
//         max
//         (
//             magSqr(boundBox(zoneB().localPoints()).span()),
//             magSqr(bbA.span())
//         )
        great
    );
    const vectorField& pCf = zoneB().localPoints();
    forAll(pCf, facei)
    {
        const vector& pt = pCf[facei];
        pointIndexHit pIH = tree.findNearest(pt, nds);
        if (pIH.hit())
        {
            zoneToZoneMap[facei] = pIH.index();
        }
    }

    if (requireMatch_ && gMin(zoneToZoneMap) == -1)
    {
        FatalErrorInFunction
            << "Cannot calculate the map between interfaces!" << nl
            << "ZoneA and zoneB are not similar (zoneA patch ="
            << " " << globalPatchA().patchName() << ", zoneB "
            << "patch = " << globalPatchB().patchName() << ")" << endl
            << abort(FatalError);
    }
}


const labelList&
nearestNeighbourPatchToPatchMapping::zoneAToZoneBPointMap() const
{
    if (zoneAToZoneBPointMapPtr_.empty())
    {
        calcZoneAToZoneBPointMap();
    }

    return zoneAToZoneBPointMapPtr_;
}


void nearestNeighbourPatchToPatchMapping::calcZoneBToZoneAPointMap() const
{
     if (zoneBToZoneAPointMapPtr_.valid())
    {
        FatalErrorInFunction
            << "List already set!" << abort(FatalError);
    }

    // Check zones are conformal
    checkZoneSizes();

    // Initialise map
    zoneBToZoneAPointMapPtr_.set
    (
        new labelList(labelList(zoneA().nPoints(), -1))
    );
    labelList& zoneToZoneMap = zoneBToZoneAPointMapPtr_();

    treeBoundBox bbB(zoneB().localPoints());
    bbB = bbB.extend(1e-4);

    indexedOctree<treeDataPoint> tree
    (
        treeDataPoint(zoneB().localPoints()),
        bbB,
        10,
        10,
        3
    );

    const scalar nds
    (
//         max
//         (
//             magSqr(boundBox(zoneA().localPoints()).span()),
//             magSqr(bbB.span())
//         )
        great
    );
    const vectorField& pCf = zoneA().localPoints();
    forAll(pCf, facei)
    {
        const vector& pt = pCf[facei];
        pointIndexHit pIH = tree.findNearest(pt, nds);
        if (pIH.hit())
        {
            zoneToZoneMap[facei] = pIH.index();
        }
    }

    if (requireMatch_ && gMin(zoneToZoneMap) == -1)
    {
        FatalErrorInFunction
            << "Cannot calculate the map between interfaces!" << nl
            << "ZoneA and zoneB are not similar (zoneA patch ="
            << " " << globalPatchA().patchName() << ", zoneB "
            << "patch = " << globalPatchB().patchName() << ")" << endl
            << abort(FatalError);
    }
}


const labelList&
nearestNeighbourPatchToPatchMapping::zoneBToZoneAPointMap() const
{
    if (zoneBToZoneAPointMapPtr_.empty())
    {
        calcZoneBToZoneAPointMap();
    }

    return zoneBToZoneAPointMapPtr_();

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nearestNeighbourPatchToPatchMapping::nearestNeighbourPatchToPatchMapping
(
    const dictionary& dict,
    const primitivePatch& patchA,
    const primitivePatch& patchB,
    const globalPolyPatch& globalPatchA,
    const globalPolyPatch& globalPatchB
)
:
    patchToPatchMapping
    (
        typeName_(), dict, patchA, patchB, globalPatchA, globalPatchB
    ),
    zoneAToZoneBFaceMapPtr_(),
    zoneBToZoneAFaceMapPtr_(),
    zoneAToZoneBPointMapPtr_(),
    zoneBToZoneAPointMapPtr_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

labelList nearestNeighbourPatchToPatchMapping::unmappedFaces
(
    const globalPolyPatch& patch
) const
{
    const labelList& addr =
        (&patch == &(globalPatchA()))
      ? zoneBToZoneAFaceMap()
      : zoneAToZoneBFaceMap();
    DynamicList<label> unmapped;
    forAll(addr, i)
    {
        if (addr[i] < 0)
        {
            unmapped.append(i);
        }
    }
    return unmapped;
}


labelList nearestNeighbourPatchToPatchMapping::unmappedPoints
(
    const globalPolyPatch& patch
) const
{
    const labelList& addr =
        (&patch == &(globalPatchA()))
      ? zoneBToZoneAPointMap()
      : zoneAToZoneBPointMap();
    DynamicList<label> unmapped;
    forAll(addr, i)
    {
        if (addr[i] < 0)
        {
            unmapped.append(i);
        }
    }
    return unmapped;
}


void nearestNeighbourPatchToPatchMapping::transferFaces
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<scalar>& fromField,  // from field
    Field<scalar>& toField           // to field
) const
{
    transferFaces<scalar>(fromZone, toZone, fromField, toField);
}


void nearestNeighbourPatchToPatchMapping::transferFaces
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<vector>& fromField,  // from field
    Field<vector>& toField           // to field
) const
{
    transferFaces<vector>(fromZone, toZone, fromField, toField);
}


void nearestNeighbourPatchToPatchMapping::transferFaces
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<symmTensor>& fromField,  // from field
    Field<symmTensor>& toField           // to field
) const
{
    transferFaces<symmTensor>(fromZone, toZone, fromField, toField);
}


void nearestNeighbourPatchToPatchMapping::transferFaces
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<sphericalTensor>& fromField,  // from field
    Field<sphericalTensor>& toField           // to field
) const
{
    transferFaces<sphericalTensor>(fromZone, toZone, fromField, toField);
}


void nearestNeighbourPatchToPatchMapping::transferFaces
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<tensor>& fromField,  // from field
    Field<tensor>& toField           // to field
) const
{
    transferFaces<tensor>(fromZone, toZone, fromField, toField);
}


void nearestNeighbourPatchToPatchMapping::transferPoints
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<scalar>& fromField,  // from field
    Field<scalar>& toField           // to field
) const
{
    transferPoints<scalar>(fromZone, toZone, fromField, toField);
}

void nearestNeighbourPatchToPatchMapping::transferPoints
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<vector>& fromField,  // from field
    Field<vector>& toField           // to field
) const
{
    transferPoints<vector>(fromZone, toZone, fromField, toField);
}


void nearestNeighbourPatchToPatchMapping::transferPoints
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<symmTensor>& fromField,  // from field
    Field<symmTensor>& toField           // to field
) const
{
    transferPoints<symmTensor>(fromZone, toZone, fromField, toField);
}


void nearestNeighbourPatchToPatchMapping::transferPoints
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<sphericalTensor>& fromField,  // from field
    Field<sphericalTensor>& toField           // to field
) const
{
    transferPoints<sphericalTensor>(fromZone, toZone, fromField, toField);
}


void nearestNeighbourPatchToPatchMapping::transferPoints
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<tensor>& fromField,  // from field
    Field<tensor>& toField           // to field
) const
{
    transferPoints<tensor>(fromZone, toZone, fromField, toField);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchToPatchMappings

} // End namespace Foam

// ************************************************************************* //
