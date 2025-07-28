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

#include "intersectionPatchToPatchMapping.H"
#include "triIntersect.H"
#include "vtkWritePolyData.H"
#include "patchToPatchTools.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace patchToPatchMappings
{
    defineTypeNameAndDebug(intersection, 0);
    addToRunTimeSelectionTable
    (
        patchToPatchMapping,
        intersection,
        dictionary
    );

    int intersection::debugSrcFacei =
        debug::debugSwitch((intersection::typeName + "SrcFace").c_str(), -1);
    int intersection::debugTgtFacei =
        debug::debugSwitch((intersection::typeName + "TgtFace").c_str(), -1);
}
}

// * * * * * * * * * * * Public Static Member Functions * * * * * * * * * * * //

void Foam::patchToPatchMappings::intersection::generatePointWeights
(
    const primitivePatch& patch,
    const primitivePatch& otherPatch,
    const List<DynamicList<label>>& faceAddr,
    List<DynamicList<label>>& pointAddr,
    List<DynamicList<scalar>>& weights
)
{

    const pointField& points = patch.localPoints();
    const pointField& otherPoints = otherPatch.localPoints();

    labelHashSet checkedFaces;
    labelHashSet addedPoints;
    weights.setSize(points.size());
    forAll(points, pointi)
    {
        weights[pointi].clear();
        pointAddr[pointi].clear();

        // Collect all relevant points
        const point& pt = points[pointi];
        checkedFaces.clear();
        addedPoints.clear();

        label nearestOther = -1;
        scalar nearestDist = great;

        // Find the face that is closest to the point
        const labelList& pointFaces = patch.pointFaces()[pointi];
        forAll(pointFaces, pfi)
        {
            const label facei = pointFaces[pfi];
            const labelList& otherFaces = faceAddr[facei];
            forAll(otherFaces, ofi)
            {
                const label otherFacei = otherFaces[ofi];
                if (checkedFaces.insert(otherFacei))
                {
                    pointHit ph = otherPatch[otherFacei].nearestPoint
                    (
                        pt,
                        otherPoints
                    );

                    // No check if hit since the faces are correctly mapped
                    // so the nearest point should be valid
                    if (ph.distance() < nearestDist)
                    {
                        nearestOther = otherFacei;
                        nearestDist = ph.distance();
                    }
                }
            }
        }

        // If there is a face overlap, set the weights using the
        // triangle face areas
        if (nearestOther >= 0)
        {
            const face& otherFace = otherPatch[nearestOther];
            forAll(otherFace, pi)
            {
                const label p0 = otherFace.rcIndex(pi);
                const label p1 = otherFace.fcIndex(pi);
                const scalar At =
                    triPointRef
                    (
                        otherPoints[otherFace[p0]],
                        otherPoints[otherFace[pi]],
                        otherPoints[otherFace[p1]]
                    ).mag();
                const scalar A0 =
                    triPointRef
                    (
                        otherPoints[otherFace[p0]],
                        otherPoints[otherFace[pi]],
                        pt
                    ).mag();
                const scalar A1 =
                    triPointRef
                    (
                        otherPoints[otherFace[pi]],
                        otherPoints[otherFace[p1]],
                        pt
                    ).mag();

                const scalar w = At/max(A0*A1, vSmall);
                pointAddr[pointi].append(otherFace[pi]);
                weights[pointi].append(w);
            }
        }
    }
}


void Foam::patchToPatchMappings::intersection::generateFaceWeights
(
    const List<DynamicList<couple>>& couples,
    const List<scalar>& coverage,
    List<DynamicList<scalar>>& weights
)
{
    weights.setSize(couples.size());
    forAll(couples, facei)
    {
        weights[facei].resize(couples[facei].size());
        scalar aSum = 0;

        forAll(couples[facei], i)
        {
            const scalar a = mag(couples[facei][i].area);
            weights[facei][i] = a;
            aSum += a;
        }

        forAll(couples[facei], i)
        {
            weights[facei][i] *= coverage[facei]/aSum;
        }
    }
}


// * * * * * * * * * * * Private Static Member Functions * * * * * * * * * * //

template<class Type>
Foam::FixedList<Type, 3>
Foam::patchToPatchMappings::intersection::triPointValues
(
    const triFace& t,
    const UList<Type>& values
)
{
    FixedList<Type, 3> result;
    forAll(t, i)
    {
        result[i] = values[t[i]];
    }
    return result;
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::treeBoundBox Foam::patchToPatchMappings::intersection::makeBb
(
    const face& f,
    const pointField& pts,
    const vectorField& pns
) const
{
    DynamicList<point> ps;
    ps.clear();

    const scalar l = sqrt(mag(f.area(pts)));
    forAll(f, fpi)
    {
        const label pointi = f[fpi];

        const point& pt = pts[pointi];
        const vector& n = pns[pointi];

        ps.append(pt - 0.5*(l*n));
        ps.append(pt + 0.5*(l*n));
    }
    return treeBoundBox(ps);
}

bool Foam::patchToPatchMappings::intersection::intersectFaces
(
    const primitivePatch& srcPatch,
    const pointField& srcPts0,
    const primitivePatch& tgtPatch,
    const pointField& tgtPts0,
    const vectorField& pointNormals,
    const vectorField& pointNormals0,
    const label srcFacei,
    const label tgtFacei
)
{
    const treeBoundBox srcFaceBb = makeBb
    (
        srcPatch.localFaces()[srcFacei],
        srcPatch.localPoints(),
        pointNormals
    );
    const treeBoundBox tgtFaceBb(tgtPatch.points(), tgtPatch[tgtFacei]);

    if (!srcFaceBb.overlaps(tgtFaceBb))
    {
        return false;
    }

    if (srcTriPoints_[srcFacei].empty())
    {
        triEngine_.triangulate
        (
            UIndirectList<point>
            (
                srcPatch.localPoints(),
                srcPatch.localFaces()[srcFacei]
            )
        );
        srcTriPoints_[srcFacei] =
            triEngine_.triPoints(srcPatch.localFaces()[srcFacei]);
        srcTriFaceEdges_[srcFacei] = triEngine_.triEdges();
    }
    if (tgtTriPoints_[tgtFacei].empty())
    {
        triEngine_.triangulate
        (
            UIndirectList<point>
            (
                tgtPatch.localPoints(),
                tgtPatch.localFaces()[tgtFacei]
            )
        );
        tgtTriPoints_[tgtFacei] =
            triEngine_.triPoints(tgtPatch.localFaces()[tgtFacei]);
        tgtTriFaceEdges_[tgtFacei] = triEngine_.triEdges();
    }

    // Construct and initialise workspace
    bool srcCouples = false;
    couple srcCouple;
    srcFaceEdgePart_.resize(srcPatch[srcFacei].size());
    forAll(srcFaceEdgePart_, srcFaceEdgei)
    {
        const edge e =
            srcPatch.localFaces()[srcFacei].faceEdge(srcFaceEdgei);
        const vector eC = e.centre(srcPatch.localPoints());
        srcFaceEdgePart_[srcFaceEdgei] = part(Zero, eC);
    }

    bool tgtCouples = false;
    couple tgtCouple;
    tgtFaceEdgePart_.resize(tgtPatch[tgtFacei].size());
    forAll(tgtFaceEdgePart_, tgtFaceEdgei)
    {
        const edge e =
            tgtPatch.localFaces()[tgtFacei].faceEdge(tgtFaceEdgei);
        const vector eC = e.centre(tgtPatch.localPoints());
        tgtFaceEdgePart_[tgtFaceEdgei] = part(Zero, eC);
    }

    part errorPart(Zero, srcPatch.faceCentres()[srcFacei]);

    // Cache the face area magnitudes
    const scalar srcMagA = mag(srcPatch.faceAreas()[srcFacei]);
    const scalar tgtMagA = mag(tgtPatch.faceAreas()[tgtFacei]);

    // Determine whether or not to debug this tri intersection
    const bool debugTriIntersect =
        (debugSrcFacei != -1 || debugTgtFacei != -1)
     && (debugSrcFacei == -1 || debugSrcFacei == srcFacei)
     && (debugTgtFacei == -1 || debugTgtFacei == tgtFacei);

    // Loop the face triangles and compute the intersections
    bool anyCouples = false;
    forAll(srcTriPoints_[srcFacei], srcFaceTrii)
    {
        const triFace& srcT = srcTriPoints_[srcFacei][srcFaceTrii];

        const FixedList<point, 3> srcPs
        (
            triPointValues(srcT, srcPatch.localPoints())
        );
        const FixedList<vector, 3> srcNs(triPointValues(srcT, pointNormals));

        forAll(tgtTriPoints_[tgtFacei], tgtFaceTrii)
        {
            const triFace tgtT =
                reverse()
              ? tgtTriPoints_[tgtFacei][tgtFaceTrii].reverseFace()
              : tgtTriPoints_[tgtFacei][tgtFaceTrii];

            const FixedList<point, 3> tgtPs =
                triPointValues(tgtT, tgtPatch.localPoints());

            // Do tri-intersection
            ictSrcPoints_.clear();
            ictSrcPointNormals_.clear();
            ictTgtPoints_.clear();
            ictPointLocations_.clear();
            triIntersect::intersectTris
            (
                srcPs,
                srcNs,
                {false, false, false},
                {-1, -1, -1},
                tgtPs,
                {false, false, false},
                {-1, -1, -1},
                ictSrcPoints_,
                ictSrcPointNormals_,
                ictTgtPoints_,
                ictPointLocations_,
                debugTriIntersect,
                debugTriIntersect
              ? word
                (
                    typeName
                  + "_srcFace=" + Foam::name(srcFacei)
                  + "_tgtFace=" + Foam::name(tgtFacei)
                  + "_intersection=" + Foam::name
                    (srcFaceTrii*tgtTriPoints_[tgtFacei].size() + tgtFaceTrii)
                )
              : word::null
            );

            // If there is no intersection then continue
            if (ictPointLocations_.empty())
            {
                continue;
            }

            // Mark that there has been an intersection
            anyCouples = true;

            // Compute the intersection geometry
            const part ictSrcPart(ictSrcPoints_);
            const part ictTgtPart(ictTgtPoints_);

            // If the intersection is below tolerance then continue
            if
            (
                mag(ictSrcPart.area) < small*srcMagA
             || mag(ictTgtPart.area) < small*tgtMagA
            )
            {
                continue;
            }

            // Mark that the source and target faces intersect
            srcCouples = tgtCouples = true;

            // Store the intersection geometry
            srcCouple += ictSrcPart;
            srcCouple.nbr += ictTgtPart;
            if (reverse())
            {
                tgtCouple += ictTgtPart;
                tgtCouple.nbr += ictSrcPart;
            }
            else
            {
                tgtCouple -= ictTgtPart;
                tgtCouple.nbr -= ictSrcPart;
            }

            // Store the intersection polygons for debugging
            const label debugSrcPoint0 = debugPoints_.size();
            const label debugTgtPoint0 =
                debugPoints_.size() + ictSrcPoints_.size();
            if (debug)
            {
                debugPoints_.append(ictSrcPoints_);
                debugPoints_.append(ictTgtPoints_);
                debugFaces_.append
                (
                    debugSrcPoint0 + identityMap(ictSrcPoints_.size())
                );
                debugFaceSrcFaces_.append(srcFacei);
                debugFaceTgtFaces_.append(tgtFacei);
                debugFaceSides_.append(1);
                debugFaces_.append
                (
                    debugTgtPoint0 + identityMap(ictTgtPoints_.size())
                );
                debugFaceSrcFaces_.append(srcFacei);
                debugFaceTgtFaces_.append(tgtFacei);
                debugFaceSides_.append(-1);
            }

            // Store edge and error areas
            forAll(ictPointLocations_, i0)
            {
                const label i1 = ictPointLocations_.fcIndex(i0);

                // Get the locations on each end of this edge of the
                // intersection polygon
                const triIntersect::location l0 = ictPointLocations_[i0];
                const triIntersect::location l1 = ictPointLocations_[i1];

                // Get the geometry for the projection of this edge
                const part ictEdgePart
                (
                    FixedList<point, 4>
                    ({
                        ictSrcPoints_[i0],
                        ictSrcPoints_[i1],
                        ictTgtPoints_[i1],
                        ictTgtPoints_[i0]
                    })
                );

                // Store the "side" of the intersection that this edge
                // corresponds to
                label ictEdgeSide = -labelMax;

                // If this edge corresponds to an edge of the source
                // triangle
                if
                (
                    l0.isSrcNotTgtPoint()
                 || l1.isSrcNotTgtPoint()
                 || (
                        l0.isIntersection()
                     && l1.isIntersection()
                     && l0.srcEdgei() == l1.srcEdgei()
                    )
                )
                {
                    const label srcEi =
                        l0.isSrcPoint() ? l0.srcPointi()
                      : l1.isSrcPoint() ? (l1.srcPointi() + 2) % 3
                      : l0.srcEdgei();

                    const label srcFaceEdgei =
                        srcTriFaceEdges_[srcFacei][srcFaceTrii][srcEi];

                    if (srcFaceEdgei < srcPatch[srcFacei].size())
                    {
                        srcFaceEdgePart_[srcFaceEdgei] += ictEdgePart;
                        ictEdgeSide = 1;
                    }
                    else
                    {
                        errorPart += ictEdgePart;
                        ictEdgeSide = 0;
                    }
                }

                // If this edge corresponds to an edge of the target
                // triangle
                else if
                (
                    l0.isTgtNotSrcPoint()
                 || l1.isTgtNotSrcPoint()
                 || (
                        l0.isIntersection()
                     && l1.isIntersection()
                     && l0.tgtEdgei() == l1.tgtEdgei()
                    )
                )
                {
                    const label tgtEi =
                        l0.isTgtPoint() ? (l0.tgtPointi() + 2) % 3
                      : l1.isTgtPoint() ? l1.tgtPointi()
                      : l0.tgtEdgei();

                    const label tgtFaceEdgei =
                        tgtTriFaceEdges_[tgtFacei][tgtFaceTrii]
                        [reverse() ? 2 - tgtEi : tgtEi];

                    if (tgtFaceEdgei < tgtPatch[tgtFacei].size())
                    {
                        if (reverse())
                        {
                            tgtFaceEdgePart_[tgtFaceEdgei] += ictEdgePart;
                        }
                        else
                        {
                            tgtFaceEdgePart_[tgtFaceEdgei] -= ictEdgePart;
                        }
                        ictEdgeSide = -1;
                    }
                    else
                    {
                        errorPart += ictEdgePart;
                        ictEdgeSide = 0;
                    }
                }

                // No other location combinations should be possible for an
                // intersection without any shared points
                else
                {
                    FatalErrorInFunction
                        << "Tri-intersection topology not recognised. "
                        << "This is a bug." << exit(FatalError);
                }

                // Store the projected edge quadrilateral for debugging
                if (debug)
                {
                    debugFaces_.append
                    (
                        labelList
                        ({
                            debugSrcPoint0 + i0,
                            debugSrcPoint0 + i1,
                            debugTgtPoint0 + i1,
                            debugTgtPoint0 + i0
                        })
                    );
                    debugFaceSrcFaces_.append(srcFacei);
                    debugFaceTgtFaces_.append(tgtFacei);
                    debugFaceSides_.append(ictEdgeSide);
                }
            }
        }
    }

    // If the source face couples the target, then store the intersection
    if (srcCouples)
    {
        localTgtFacesToSrc_[srcFacei].append(tgtFacei);
        srcCouples_[srcFacei].append(srcCouple);
    }

    // If any intersection has occurred then store the edge and error parts
    if (anyCouples)
    {
        forAll(srcFaceEdgeParts_[srcFacei], srcFaceEdgei)
        {
            srcFaceEdgeParts_[srcFacei][srcFaceEdgei] +=
                srcFaceEdgePart_[srcFaceEdgei];
        }
        srcErrorParts_[srcFacei] -= sum(tgtFaceEdgePart_);
        srcErrorParts_[srcFacei] += errorPart;
    }

    // If the target face couples the source, then store in the intersection
    if (tgtCouples)
    {
        localSrcFacesToTgt_[tgtFacei].append(srcFacei);
        tgtCouples_[tgtFacei].append(tgtCouple);
    }

    return anyCouples;
}


void Foam::patchToPatchMappings::intersection::initialise
(
    const primitivePatch& srcPatch,
    const pointField& srcPts0,
    const primitivePatch& tgtPatch,
    const pointField& tgtPts0,
    const vectorField& pointNormals,
    const vectorField& pointNormals0
)
{
    patchToPatchMapping::initialise
    (
        srcPatch,
        srcPts0,
        tgtPatch,
        tgtPts0,
        pointNormals,
        pointNormals0
    );

    srcCouples_.resize(srcPatch.size());
    forAll(localTgtFacesToSrc_, i)
    {
        srcCouples_[i].clear();
    }

    srcEdgeParts_.resize(srcPatch.nEdges());
    forAll(srcEdgeParts_, srcEdgei)
    {
        const edge& e = srcPatch.edges()[srcEdgei];
        const point c = e.centre(srcPatch.localPoints());
        srcEdgeParts_[srcEdgei] = part(Zero, c);
    }

    srcErrorParts_.resize(srcPatch.size());
    forAll(srcErrorParts_, srcFacei)
    {
        srcErrorParts_[srcFacei] =
            part(Zero, srcPatch.faceCentres()[srcFacei]);
    }

    tgtCouples_.resize(tgtPatch.size());
    forAll(localSrcFacesToTgt_, i)
    {
        tgtCouples_[i].clear();
    }

    srcTriPoints_ = List<triFaceList>(srcPatch.size());
    srcTriFaceEdges_ = List<List<FixedList<label, 3>>>(srcPatch.size());
    tgtTriPoints_ = List<triFaceList>(tgtPatch.size());
    tgtTriFaceEdges_ = List<List<FixedList<label, 3>>>(tgtPatch.size());

    srcFaceEdgeParts_.resize(srcPatch.size());
    forAll(srcFaceEdgeParts_, srcFacei)
    {
        srcFaceEdgeParts_[srcFacei].resize(srcPatch[srcFacei].size());
        forAll(srcFaceEdgeParts_[srcFacei], srcFaceEdgei)
        {
            const label srcEdgei =
                srcPatch.faceEdges()[srcFacei][srcFaceEdgei];
            srcFaceEdgeParts_[srcFacei][srcFaceEdgei] =
                srcEdgeParts_[srcEdgei];
        }
    }

    if (debug)
    {
        debugPoints_.clear();
        debugFaces_.clear();
        debugFaceSrcFaces_.clear();
        debugFaceTgtFaces_.clear();
        debugFaceSides_.clear();
    }
}


void Foam::patchToPatchMappings::intersection::generatePointWeights
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch
)
{
    generatePointWeights
    (
        srcPatch,
        tgtPatch,
        localTgtFacesToSrc_,
        localTgtPointsToSrc_,
        srcPointWeights_
    );
    generatePointWeights
    (
        tgtPatch,
        srcPatch,
        localSrcFacesToTgt_,
        localSrcPointsToTgt_,
        tgtPointWeights_
    );
}


void Foam::patchToPatchMappings::intersection::generateFaceWeights()
{
    generateFaceWeights
    (
        srcCouples_,
        srcCoverage_,
        srcFaceWeights_
    );
    generateFaceWeights
    (
        tgtCouples_,
        tgtCoverage_,
        tgtFaceWeights_
    );
}


Foam::labelList Foam::patchToPatchMappings::intersection::finaliseLocalPoints
(
    const primitivePatch& srcPatch,
    const pointField& srcPts0,
    const primitivePatch& tgtPatch,
    const pointField& tgtPts0,
    const vectorField& pointNormals,
    const vectorField& pointNormals0
)
{
    generatePointWeights(srcPatch, tgtPatch);

    const labelList newToOld
    (
        patchToPatchMapping::finaliseLocalPoints
        (
            srcPatch,
            srcPts0,
            tgtPatch,
            tgtPts0,
            pointNormals,
            pointNormals0
        )
    );
    tgtPointWeights_ = List<DynamicList<scalar>>(tgtPointWeights_, newToOld);

    return newToOld;
}


Foam::labelList Foam::patchToPatchMappings::intersection::finaliseLocalFaces
(
    const primitivePatch& srcPatch,
    const pointField& srcPts0,
    const primitivePatch& tgtPatch,
    const pointField& tgtPts0,
    const vectorField& pointNormals,
    const vectorField& pointNormals0
)
{
    const labelList newToOld =
        patchToPatchMapping::finaliseLocalFaces
        (
            srcPatch,
            srcPts0,
            tgtPatch,
            tgtPts0,
            pointNormals,
            pointNormals0
        );

    tgtCouples_ = List<DynamicList<couple>>(tgtCouples_, newToOld);

    return newToOld;
}


void Foam::patchToPatchMappings::intersection::rDistributeTgt
(
    const primitivePatch& tgtPatch,
    const pointField& tgtPts0
)
{
    patchToPatchMapping::rDistributeTgt(tgtPatch, tgtPts0);

    patchToPatchTools::rDistributeListList
    (
        tgtPatch.size(),
        tgtFacesMapPtr_(),
        tgtCouples_
    );

    if (needPoints_)
    {
        // Reverse distribute the point weights
        patchToPatchTools::rDistributeListList
        (
            tgtPatch.nPoints(),
            tgtPointsMapPtr_(),
            tgtPointWeights_
        );

        // Remove zero weights
        DynamicList<label> newIs;
        DynamicList<scalar> newWs;
        forAll(srcPointWeights_, srcPointi)
        {
            newIs.clear();
            newWs.clear();
            labelList& Is = localTgtPointsToSrc_[srcPointi];
            scalarList& ws = srcPointWeights_[srcPointi];
            forAll(ws, i)
            {
                if (ws[i] > vSmall)
                {
                    newIs.append(Is[i]);
                    newWs.append(ws[i]);
                }
            }
            Is = newIs;
            ws = newWs;
        }
        forAll(tgtPointWeights_, tgtPointi)
        {
            newIs.clear();
            newWs.clear();
            labelList& Is = localSrcPointsToTgt_[tgtPointi];
            scalarList& ws = tgtPointWeights_[tgtPointi];
            forAll(ws, i)
            {
                if (ws[i] > vSmall)
                {
                    newIs.append(Is[i]);
                    newWs.append(ws[i]);
                }
            }
            Is = newIs;
            ws = newWs;
        }
    }
}


Foam::label Foam::patchToPatchMappings::intersection::finalisePoints
(
    const primitivePatch& srcPatch,
    const pointField& srcPts0,
    const primitivePatch& tgtPatch,
    const pointField& tgtPts0,
    const vectorField& pointNormals,
    const vectorField& pointNormals0,
    const transformer& tgtToSrc
)
{
    if (isSingleProcess())
    {
        generatePointWeights(srcPatch, tgtPatch);
    }

    const label nCouples =
        patchToPatchMapping::finalisePoints
        (
            srcPatch,
            srcPts0,
            tgtPatch,
            tgtPts0,
            pointNormals,
            pointNormals0,
            tgtToSrc
        );

    // Normalize weights
    forAll(srcPointWeights_, srcPointi)
    {
        scalarList& ws = srcPointWeights_[srcPointi];
        const scalar sumW = max(sum(ws), vSmall);
        forAll(ws, i)
        {
            ws[i] /= sumW;
        }
    }
    forAll(tgtPointWeights_, tgtPointi)
    {
        scalarList& ws = tgtPointWeights_[tgtPointi];
        const scalar sumW = max(sum(ws), vSmall);
        forAll(ws, i)
        {
            ws[i] /= sumW;
        }
    }

    return nCouples;
}


Foam::label Foam::patchToPatchMappings::intersection::finaliseFaces
(
    const primitivePatch& srcPatch,
    const pointField& srcPts0,
    const primitivePatch& tgtPatch,
    const pointField& tgtPts0,
    const vectorField& pointNormals,
    const vectorField& pointNormals0,
    const transformer& tgtToSrc
)
{
    const label nCouples =
        patchToPatchMapping::finaliseFaces
        (
            srcPatch,
            srcPts0,
            tgtPatch,
            tgtPts0,
            pointNormals,
            pointNormals0,
            tgtToSrc
        );

    // Convert face-edge-parts to edge-parts
    labelList srcEdgeNParts(srcEdgeParts_.size(), 0);
    forAll(srcEdgeParts_, srcEdgei)
    {
        const edge& e = srcPatch.edges()[srcEdgei];

        srcEdgeParts_[srcEdgei] = part();

        forAll(srcPatch.edgeFaces()[srcEdgei], i)
        {
            const label srcFacei = srcPatch.edgeFaces()[srcEdgei][i];
            const label srcFaceEdgei =
                findIndex(srcPatch.faceEdges()[srcFacei], srcEdgei);

            const edge fe =
                srcPatch.localFaces()[srcFacei].faceEdge(srcFaceEdgei);

            if (edge::compare(e, fe) > 0)
            {
                srcEdgeParts_[srcEdgei] +=
                    srcFaceEdgeParts_[srcFacei][srcFaceEdgei];
            }
            else
            {
                srcEdgeParts_[srcEdgei] -=
                    srcFaceEdgeParts_[srcFacei][srcFaceEdgei];
            }

            srcEdgeNParts[srcEdgei] ++;
        }
    }
    forAll(srcEdgeParts_, srcEdgei)
    {
        srcEdgeParts_[srcEdgei].area /= srcEdgeNParts[srcEdgei];
    }

    // Add the difference between the face-edge-part and the edge-part into the
    // face-error-parts
    forAll(srcEdgeParts_, srcEdgei)
    {
        const edge& e = srcPatch.edges()[srcEdgei];

        forAll(srcPatch.edgeFaces()[srcEdgei], i)
        {
            const label srcFacei = srcPatch.edgeFaces()[srcEdgei][i];
            const label srcFaceEdgei =
                findIndex(srcPatch.faceEdges()[srcFacei], srcEdgei);

            const edge fe =
                srcPatch.localFaces()[srcFacei].faceEdge(srcFaceEdgei);

            if (edge::compare(e, fe) > 0)
            {
                srcErrorParts_[srcFacei] -= srcEdgeParts_[srcEdgei];
            }
            else
            {
                srcErrorParts_[srcFacei] += srcEdgeParts_[srcEdgei];
            }

            srcErrorParts_[srcFacei] +=
                srcFaceEdgeParts_[srcFacei][srcFaceEdgei];
        }
    }

    // Transform the target couples back to the target side
    if (!isNull(tgtToSrc))
    {
        forAll(tgtCouples_, tgtFacei)
        {
            forAll(tgtCouples_[tgtFacei], i)
            {
                couple& c = tgtCouples_[tgtFacei][i];

                c.area = tgtToSrc.invTransform(c.area);
                c.centre = tgtToSrc.invTransformPosition(c.centre);
                c.nbr.area = tgtToSrc.invTransform(c.nbr.area);
                c.nbr.centre = tgtToSrc.invTransformPosition(c.nbr.centre);
            }
        }
    }

    // Calculate coverage and total areas on both sides
    auto coverage = []
    (
        const primitivePatch& patch,
        const List<DynamicList<couple>>& couples,
        scalar& area,
        scalar& coupleArea,
        List<scalar>& coverage
    )
    {
        area = 0;
        coupleArea = 0;
        coverage.resize(patch.size());

        forAll(patch, facei)
        {
            const scalar magA = mag(patch.faceAreas()[facei]);

            vector aCouple = Zero;
            forAll(couples[facei], i)
            {
                aCouple += couples[facei][i].area;
            }
            const scalar magACouple = mag(aCouple);

            area += magA;
            coupleArea += magACouple;
            coverage[facei] = magACouple/magA;
        }

        reduce(area, sumOp<scalar>());
        reduce(coupleArea, sumOp<scalar>());
    };
    scalar srcArea = 0, srcCoupleArea = 0;
    scalar tgtArea = 0, tgtCoupleArea = 0;
    coverage(srcPatch, srcCouples_, srcArea, srcCoupleArea, srcCoverage_);
    coverage(tgtPatch, tgtCouples_, tgtArea, tgtCoupleArea, tgtCoverage_);

    // Clear the triangulation workspace
    srcTriPoints_.clear();
    srcTriFaceEdges_.clear();
    tgtTriPoints_.clear();
    tgtTriFaceEdges_.clear();

    // Clear face-edge-parts
    srcFaceEdgePart_.clear();
    tgtFaceEdgePart_.clear();
    srcFaceEdgeParts_.clear();

    // Checking and reporting
    if (nCouples != 0)
    {
        scalarField srcOpenness(srcPatch.size());
        scalarField srcError(srcPatch.size());
        scalarField srcDepth(srcPatch.size());
        scalarField srcAngle(srcPatch.size());
        forAll(srcPatch, srcFacei)
        {
            const vector& a = srcPatch.faceAreas()[srcFacei];
            const scalar magA = mag(a);
            const point& c = srcPatch.faceCentres()[srcFacei];

            couple Cpl(part(Zero, c), part(Zero, c));
            forAll(srcCouples_[srcFacei], srcTgtFacei)
            {
                const couple& cpl = srcCouples_[srcFacei][srcTgtFacei];

                Cpl += cpl;
                Cpl.nbr += cpl.nbr;
            }

            vector projectionA = Zero;
            scalar projectionV = 0;
            forAll(srcCouples_[srcFacei], srcTgtFacei)
            {
                const couple& cpl = srcCouples_[srcFacei][srcTgtFacei];

                projectionA += cpl.nbr.area;
                projectionV +=
                    - (cpl.area/3 & (cpl.centre - Cpl.centre))
                    + (cpl.nbr.area/3 & (cpl.nbr.centre - Cpl.centre));
            }
            forAll(srcPatch.faceEdges()[srcFacei], srcFaceEdgei)
            {
                const label srcEdgei =
                    srcPatch.faceEdges()[srcFacei][srcFaceEdgei];

                const edge& e = srcPatch.edges()[srcEdgei];
                const edge fe =
                    srcPatch.localFaces()[srcFacei].faceEdge(srcFaceEdgei);

                const scalar sign = edge::compare(e, fe);

                projectionA += sign*srcEdgeParts_[srcEdgei].area;
                projectionV +=
                    sign*srcEdgeParts_[srcEdgei].area/3
                  & (srcEdgeParts_[srcEdgei].centre - Cpl.centre);
            }
            projectionA += srcErrorParts_[srcFacei].area;
            projectionV +=
                srcErrorParts_[srcFacei].area/3
              & (srcErrorParts_[srcFacei].centre - Cpl.centre);

            const vector aHat = normalised(a);
            const vector aOppHat = normalised(a - Cpl.area + Cpl.nbr.area);
            srcAngle[srcFacei] =
                radToDeg(acos(min(max(aHat & aOppHat, -1), +1)));
            srcOpenness[srcFacei] = mag(projectionA - Cpl.area)/magA;
            srcError[srcFacei] = mag(srcErrorParts_[srcFacei].area)/magA;
            srcDepth[srcFacei] = mag(projectionV)/pow3(sqrt(magA));
        }

        reduce(tgtArea, sumOp<scalar>());
        reduce(tgtCoupleArea, sumOp<scalar>());

        Info<< indent << "Source min/average/max coverage = "
            << gMin(srcCoverage_) << '/' << srcCoupleArea/srcArea << '/'
            << gMax(srcCoverage_) << endl
            << indent << "Target min/average/max coverage = "
            << gMin(tgtCoverage_) << '/' << tgtCoupleArea/tgtArea << '/'
            << gMax(tgtCoverage_) << endl
            << indent << "Source average openness/error/depth/angle = "
            << gAverage(srcOpenness) << '/' << gAverage(srcError) << '/'
            << gAverage(srcDepth) << '/' << gAverage(srcAngle) << endl
            << indent << "Source max openness/error/depth/angle = "
            << gMax(srcOpenness) << '/' << gMax(srcError) << '/'
            << gMax(srcDepth) << '/' << gMax(srcAngle) << endl;

        if (debug)
        {
            word name = patchToPatchMapping::typeName + '_' + typeName;

            if (Pstream::parRun())
            {
                name += "_proc" + Foam::name(Pstream::myProcNo());
            }

            Info<< indent << "Writing intersected faces to "
                << name + ".vtk" << endl;
            vtkWritePolyData::write
            (
                name + ".vtk",
                name,
                false,
                debugPoints_,
                labelList(),
                labelListList(),
                debugFaces_,
                "srcFace", false, Field<label>(debugFaceSrcFaces_),
                "tgtFace", false, Field<label>(debugFaceTgtFaces_),
                "side", false, Field<label>(debugFaceSides_)
            );

            debugPoints_.clear();
            debugFaces_.clear();
            debugFaceSrcFaces_.clear();
            debugFaceTgtFaces_.clear();
            debugFaceSides_.clear();

            Info<< indent << "Writing source patch to "
                << name + "_srcPatch.vtk" << endl;
            vtkWritePolyData::write
            (
                name + "_srcPatch" + ".vtk",
                name + "_srcPatch",
                false,
                srcPatch.localPoints(),
                labelList(),
                labelListList(),
                srcPatch.localFaces(),
                "coverage", false, scalarField(srcCoverage_),
                "openness", false, srcOpenness,
                "error", false, srcError,
                "depth", false, srcDepth,
                "angle", false, srcAngle,
                "normals", true, pointNormals
            );

            Info<< indent << "Writing target patch to "
                << name + "_tgtPatch.vtk" << endl;
            vtkWritePolyData::write
            (
                name + "_tgtPatch" + ".vtk",
                name + "_tgtPatch",
                false,
                tgtPatch.localPoints(),
                labelList(),
                labelListList(),
                tgtPatch.localFaces(),
                "coverage", false, scalarField(tgtCoverage_)
            );
        }
    }

    generateFaceWeights();

    return nCouples;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchToPatchMappings::intersection::intersection
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch,
    const dictionary& dict,
    const bool needPoints,
    const bool reverse
)
:
    patchToPatchMapping(srcPatch, tgtPatch, dict, needPoints, reverse)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
