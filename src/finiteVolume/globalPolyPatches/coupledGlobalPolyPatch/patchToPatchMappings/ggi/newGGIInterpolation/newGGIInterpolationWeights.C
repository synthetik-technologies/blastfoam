/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

Description
    Mass-conservative face interpolation of face data between two
    primitivePatches

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

Modification by:
    Martin Beaudoin, Hydro-Quebec, (2008)

\*---------------------------------------------------------------------------*/

#include "newGGIInterpolationTemplate.H"
#include "objectHit.H"
#include "boolList.H"
#include "DynamicList.H"
#include "dimensionedConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
const Foam::debug::tolerancesSwitch
newGGIInterpolation<MasterPatch, SlavePatch>::areaErrorTol_
(
    "GGIAreaErrorTol",
    1.0e-8,
    "Minimum GGI face to face intersection area.  "
    "The smallest accepted GGI weighting factor."
);


template<class MasterPatch, class SlavePatch>
const Foam::debug::tolerancesSwitch
newGGIInterpolation<MasterPatch, SlavePatch>::featureCosTol_
(
    "GGIFeatureCosTol",
    0.8,
    "Minimum cosine value between 2 GGI patch neighbouring facet normals."
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
void newGGIInterpolation<MasterPatch, SlavePatch>::calcAddressing() const
{
    if
    (
        masterAddrPtr_
     || masterWeightsPtr_
     || slaveAddrPtr_
     || slaveWeightsPtr_
     || uncoveredMasterAddrPtr_
     || uncoveredSlaveAddrPtr_
    )
    {
        FatalErrorIn
        (
            "void newGGIInterpolation<MasterPatch, SlavePatch>::"
            "calcAddressing() const"
        )   << "Addressing already calculated"
            << abort(FatalError);
    }

    if (debug)
    {
        InfoIn
        (
            "void newGGIInterpolation<MasterPatch, SlavePatch>::"
            "calcAddressing() const"
        )   << "Evaluation of GGI weighting factors:" << endl;
    }

    // Create the dynamic lists to hold the addressing

    // First, find a rough estimate of each slave and master facet
    // neighborhood by filtering out all the faces located outside of
    // an Axis-Aligned Bounding Box (AABB).  Warning: This algorithm
    // is based on the evaluation of AABB boxes, which is pretty fast;
    // but still the complexity of the algorithm is n^2, wich is
    // pretty bad for GGI patches composed of 100,000 of facets...  So
    // here is the place where we could certainly gain major speedup
    // for larger meshes.

    // The candidates master neighbours
    // Choice of algorithm:
    // 1) Axis-aligned bounding box
    // 2) Octree search with bounding box
    // 3) 3-D vector distance


    // Note: Allocated to local size for parallel search.  HJ, 27/Apr/2016
    labelListList candidateMasterNeighbors;

    if (usePrevCandidateMasterNeighbors_)
    {
        updateNeighboursAABB(candidateMasterNeighbors);
    }
    else if (reject_ == AABB)
    {
         findNeighboursAABB(candidateMasterNeighbors);
    }
    else if (reject_ == BB_OCTREE)
    {
         findNeighboursBBOctree(candidateMasterNeighbors);
    }
    else if (reject_ == THREE_D_DISTANCE)
    {
         findNeighbours3D(candidateMasterNeighbors);
    }
    else
    {
        FatalErrorIn
        (
            "void newGGIInterpolation<MasterPatch, SlavePatch>::"
            "calcAddressing() const"
        )   << "Unknown search"
            << abort(FatalError);
    }

    // Next, we move to the 2D world.  We project each slave and
    // master face onto a local plane defined by the master face
    // normal.  We filter out a few false neighbors using the
    // Separating Axes Theorem

    // It is in this local plane that we will refine our list of
    // neighbors.  So for a given a neighbor face, we need as many
    // projections as there are neighbors closeby.

    const pointField& masterPatchPoints = masterPatch_.localPoints();
    const vectorField masterPatchNormals = masterPatch_.faceNormals();

    // ZT, 05/07/2014
    const vectorField& slavePatchNormals = slavePatch_.faceNormals();

    // Store the polygon made by projecting the face points onto the
    // face normal
    // The master faces polygons
    List<pointField> masterFace2DPolygon(masterPatch_.size());

    // Tolerance factor for the Separation of Axes Theorem == distErrorTol_

    // The final master/slave list, after filtering out the "false" neighbours
    // Note: rescale only the local ones. HJ, 27/Apr/2016
    // Note: slave neighbours and weights delayed until master cutting
    // is complete.  HJ, 27/Apr/2016
    List<DynamicList<label> > masterNeighbors(masterPatch_.size());
    List<DynamicList<scalar> > masterNeighborsWeights(masterPatch_.size());

    // Store slave negighbour weights at the same time, but under
    // addressing.  The list will be redistributed to slaves in a separate
    // loop.  HJ, 27/Apr/2016
    List<DynamicList<scalar> > slaveOnMasterNeighborsWeights
    (
        masterPatch_.size()
    );

    // Parallel search split.  HJ, 27/Apr/2016
    const label pmStart = this->parMasterStart();

    for
    (
        label faceMi = this->parMasterStart();
        faceMi < this->parMasterEnd();
        faceMi++
    )
//     forAll (masterPatch_, faceMi)
    {
        // Set capacity
        masterNeighbors[faceMi].setCapacity(8);

        // First, we make sure that all the master faces points are
        // recomputed onto the 2D plane defined by the master faces
        // normals.
        // For triangles, this is useless, but for N-gons
        // with more than 3 points, this is essential.
        // The intersection between the master and slave faces will be
        // done in these 2D reference frames

        // A few basic information to keep close-by
        vector currentMasterFaceNormal = masterPatchNormals[faceMi];
        vector currentMasterFaceCentre =
            masterPatch_[faceMi].centre(masterPatchPoints);

        scalarField facePolygonErrorProjection;

        // Project the master faces points onto the normal face plane to
        // form a flattened polygon
        masterFace2DPolygon[faceMi] =
            projectPointsOnPlane
            (
                masterPatch_[faceMi].points(masterPatchPoints),
                currentMasterFaceCentre,
                currentMasterFaceNormal,
                facePolygonErrorProjection
            );

        // Next we compute an orthonormal basis (u, v, w) aligned with
        // the face normal for doing the 3D to 2D projection.
        //
        // "w" is aligned on the face normal.  We need to select a "u"
        // direction, it can be anything as long as it lays on the
        // projection plane.  We chose to use the direction from the
        // master face center to the most distant projected master face
        // point on the plane.  Finally, we get "v" by evaluating the
        // cross-product w^u = v.  And we make sure that u, v, and w are
        // normalized.
        //
        //
        // u  =  vector from face center to most distant projected master face point.
        //                                    /       .
        //           ^y                     / |       .      .w = normal to master face
        //           |                    /   |       .    .
        //           |                  /     |       .  .
        //           |                 |      |       .
        //           |                 |      /        .
        //           |                 |    /           .
        //           |                 |  /              .
        //           ---------> x      |/                 .
        //          /                                                 v = w^u
        //         /
        //        /
        //       z
        //
        //

        orthoNormalBasis uvw =
            computeOrthonormalBasis
            (
                currentMasterFaceCentre,
                currentMasterFaceNormal,
                masterFace2DPolygon[faceMi]
            );

        // Recompute the master polygon into this orthoNormalBasis
        // We should only see a rotation along the normal of the face here
        List<point2D> masterPointsInUV;
        scalarField masterErrorProjectionAlongW;

        masterPointsInUV =
            projectPoints3Dto2D
            (
                uvw,
                currentMasterFaceCentre,
                masterFace2DPolygon[faceMi],
                masterErrorProjectionAlongW   // Should be at zero all the way
            );

        // Compute the surface area of the polygon;
        // We need this for computing the weighting factors
        scalar surfaceAreaMasterPointsInUV = area2D(masterPointsInUV);

        // Check if polygon is CW. Should not, it should be CCW; but
        // better and cheaper to check here
        if (surfaceAreaMasterPointsInUV < 0)
        {
            reverse(masterPointsInUV);
            surfaceAreaMasterPointsInUV = -surfaceAreaMasterPointsInUV;

            // Just generate a warning until we can verify this is a non issue
            InfoIn
            (
                "void newGGIInterpolation<MasterPatch, SlavePatch>::"
                "calcAddressing()"
            )   << "The master projected polygon was CW instead of CCW.  "
                << "This is strange..."  << endl;
        }

        // Next, project the candidate master neighbours faces points
        // onto the same plane using the new orthonormal basis
        // Note: Allocated to local size for parallel search.  HJ, 27/Apr/2016
        const labelList& curCMN = candidateMasterNeighbors[faceMi - pmStart];

        forAll (curCMN, neighbI)
        {
            // For each points, compute the dot product with u,v,w.  The
            // [u,v] component will gives us the 2D cordinates we are
            // looking for for doing the 2D intersection The w component
            // is basically the projection error normal to the projection
            // plane

            // NB: this polygon is most certainly CW w/r to the uvw
            // axis because of the way the normals are oriented on
            // each side of the GGI interface... We will switch the
            // polygon to CCW in due time...
            List<point2D> neighbPointsInUV;
            scalarField neighbErrorProjectionAlongW;

            // We use the xyz points directly, with a possible transformation
            pointField curSlaveFacePoints =
                slavePatch_[curCMN[neighbI]].points(slavePatch_.localPoints());

            if (doTransform())
            {
                // Transform points to master plane
                if (forwardT_.size() == 1)
                {
                    transform
                    (
                        curSlaveFacePoints,
                        forwardT_[0],
                        curSlaveFacePoints
                    );
                }
                else
                {
                    transform
                    (
                        curSlaveFacePoints,
                        forwardT_[curCMN[neighbI]],
                        curSlaveFacePoints
                    );
                }
            }

            // Apply the translation offset in order to keep the
            // neighbErrorProjectionAlongW values to a minimum
            if (doSeparation())
            {
                if (forwardSep_.size() == 1)
                {
                    curSlaveFacePoints += forwardSep_[0];
                }
                else
                {
                    curSlaveFacePoints += forwardSep_[curCMN[neighbI]];
                }
            }

            neighbPointsInUV =
                projectPoints3Dto2D
                (
                    uvw,
                    currentMasterFaceCentre,
                    curSlaveFacePoints,
                    neighbErrorProjectionAlongW
                );

            // ZT, 05/07/2014
            // PC, re-enable 26/04/2017
            const scalar orientation =
                (
                    masterPatchNormals[faceMi]
                  & slavePatchNormals[curCMN[neighbI]]
                );

            // We are now ready to filter out the "bad" neighbours.
            // For this, we will apply the Separating Axes Theorem
            // http://en.wikipedia.org/wiki/Separating_axis_theorem.

            // This will be the second and last quick reject test.
            // We will use the 2D projected points for both the master
            // patch and its neighbour candidates
            if
            (
                detect2dPolygonsOverlap
                (
                    masterPointsInUV,
                    neighbPointsInUV,
                    sqrt(areaErrorTol_()) // distErrorTol
                )
             && (orientation < -SMALL) // ZT, 05/07/2014
            )
            {
                // We have an overlap between the master face and this
                // neighbor face.
                label faceMaster = faceMi;
                label faceSlave  = curCMN[neighbI];

                // Compute the surface area of the neighbour polygon;
                // We need this for computing the weighting factors
                scalar surfaceAreaNeighbPointsInUV = area2D(neighbPointsInUV);

                // Check for CW polygons. It most certainly is, and
                // the polygon intersection algorithms are expecting
                // to work with CCW point ordering for the polygons
                if (surfaceAreaNeighbPointsInUV < 0.0)
                {
                    reverse(neighbPointsInUV);
                    surfaceAreaNeighbPointsInUV = -surfaceAreaNeighbPointsInUV;
                }


                // We compute the intersection area using the
                // Sutherland-Hodgman algorithm.  Of course, if the
                // intersection area is 0, that would constitute the last and
                // final reject test, but it would also be an indication that
                // our 2 previous rejection tests are a bit laxed...  or that
                // maybe we are in presence of concave polygons....
                scalar intersectionArea =
                    polygonIntersection
                    (
                        masterPointsInUV,
                        neighbPointsInUV
                    );

                if (intersectionArea > VSMALL) // Or > areaErrorTol_ ???
                {
                    // We compute the GGI weights based on this
                    // intersection area, and on the individual face
                    // area on each side of the GGI.

                    // Since all the intersection have been computed
                    // in the projected UV space we need to compute
                    // the weights using the surface area from the
                    // faces projection as well. That way, we make
                    // sure all our factors will sum up to 1.0.

                    // Add slave to master
                    masterNeighbors[faceMaster].append(faceSlave);

                    // Add master weight to master
                    masterNeighborsWeights[faceMaster].append
                    (
                        intersectionArea/surfaceAreaMasterPointsInUV
                    );

                    // Record slave weight on master to avoid recalculation
                    // of projected areas.  HJ, 27/Apr/2016
                    slaveOnMasterNeighborsWeights[faceMaster].append
                    (
                        intersectionArea/surfaceAreaNeighbPointsInUV
                    );

                    // Note: Slave side will be reconstructed after the
                    // parallel cutting and reduce operations.
                    // HJ, 27/Apr/2016
                }
                else
                {
                    WarningIn
                    (
                        "newGGIInterpolation<MasterPatch, SlavePatch>::"
                        "calcAddressing()"
                    )   << "polygonIntersection is returning a "
                        << "zero surface area between " << nl
                        << "     Master face: " << faceMi
                        << " and Neighbour face: " << curCMN[neighbI]
                        << " intersection area = " << intersectionArea << nl
                        << "Please check the two quick-check algorithms for "
                        << "newGGIInterpolation.  Something is  missing."
                        << endl;
                }
            }
        }

        // We went through all the possible neighbors for this face.
    }


    // Allocate the member attributes and pack addressing
    masterAddrPtr_ = new labelListList(masterPatch_.size());
    labelListList& ma  = *masterAddrPtr_;

    // Parallel search split.  HJ, 27/Apr/2016
    for (label mfI = this->parMasterStart(); mfI < this->parMasterEnd(); mfI++)
//     forAll (ma, mfI)
    {
        ma[mfI].transfer(masterNeighbors[mfI].shrink());
    }

    // Parallel communication: reduce master addressing
    if (globalData())
    {
        Pstream::combineGather(ma, Pstream::listEq());
        Pstream::combineScatter(ma);
    }

    masterWeightsPtr_ = new scalarListList(masterPatch_.size());
    scalarListList& maW = *masterWeightsPtr_;

    // Parallel search split.  HJ, 27/Apr/2016
    for (label mfI = this->parMasterStart(); mfI < this->parMasterEnd(); mfI++)
//     forAll (maW, mfI)
    {
        maW[mfI].transfer(masterNeighborsWeights[mfI].shrink());
    }

    // Parallel communication: reduce master weights
    if (globalData())
    {
        Pstream::combineGather(maW, Pstream::listEq());
        Pstream::combineScatter(maW);
    }

    // Reduce slave on master weights
    scalarListList smaW(masterPatch_.size());

    for (label mfI = this->parMasterStart(); mfI < this->parMasterEnd(); mfI++)
//     forAll (smW, mfI)
    {
        smaW[mfI].transfer(slaveOnMasterNeighborsWeights[mfI].shrink());
    }

    // Parallel communication: reduce master weights
    if (globalData())
    {
        Pstream::combineGather(smaW, Pstream::listEq());
        Pstream::combineScatter(smaW);
    }
    // Slave neighbours and weights
    List<DynamicList<label, 8> > slaveNeighbors(slavePatch_.size());
    List<DynamicList<scalar, 8> > slaveNeighborsWeights(slavePatch_.size());

    // Note: slave side is not parallelised: book-keeping only
    // HJ, 27/Apr/2016

    // Loop through the complete patch and distribute slave neighbours
    // and weights based onthe known intersection area from the master.
    // The code has been reorgnised to allow masters cutting to be
    // performed in parallel and collect the slave side once the parallel
    // reduction is complete.  HJ, 27/Apr/2016

    // Note: loop through the complete slave side
    forAll (ma, mfI)
    {
        // Gte master neighbours and weights
        const labelList& curMa = ma[mfI];
        const scalarList& cursmaW = smaW[mfI];

        // Current master face indes
        const label faceMaster = mfI;

        forAll (curMa, mAI)
        {
            const label faceSlave = curMa[mAI];
            // Add this master as a neighbour to its slave
            slaveNeighbors[faceSlave].append(faceMaster);

            slaveNeighborsWeights[faceSlave].append(cursmaW[mAI]);
        }
    }


    slaveAddrPtr_ = new labelListList(slavePatch_.size());
    labelListList& sa = *slaveAddrPtr_;

    forAll (sa, sfI)
    {
        sa[sfI].transfer(slaveNeighbors[sfI].shrink());
    }

    slaveWeightsPtr_ = new scalarListList(slavePatch_.size());
    scalarListList& saW = *slaveWeightsPtr_;

    forAll (sa, sfI)
    {
        saW[sfI].transfer(slaveNeighborsWeights[sfI].shrink());
    }

    // Now that the neighbourhood is known, let's go hunting for
    // non-overlapping faces

    uncoveredMasterAddrPtr_ =
        new labelList
        (
            findNonOverlappingFaces(maW, masterNonOverlapFaceTol_)
        );

    // Not parallelised.  HJ, 27/Apr/2016
    uncoveredSlaveAddrPtr_ =
        new labelList
        (
            findNonOverlappingFaces(saW, slaveNonOverlapFaceTol_)
        );

    // Rescaling the weighting factors so they will sum up to 1.0
    // See the comment for the method ::rescaleWeightingFactors() for
    // more information.  By default, we always rescale.  But for some
    // special kind of GGI interpolation, like the mixingPlaneGGI,
    // then we need the brute values, so no rescaling in that
    // case. Hence the little flag rescaleGGIWeightingFactors_

    // Not parallelised.  HJ, 27/Apr/2016
    if (rescaleGGIWeightingFactors_)
    {
        rescaleWeightingFactors();
    }
}


// Rescaling the weighting factors so they will sum up to 1.0 This is
// necessary for the slave weighting factors because intersection with
// master neighbours are usually computed from different projection
// plane, so the weighting factor don't quite sum up to 1.0 For the
// slave weighting factors, we are usually talking of delta of the
// order of 10e-6 here.
//
// For the master weighting factor, this is another story. We truly
// expect that the master weighting factors will exactly sum up to 1.0
// if all the neighbours are properly identified.
//
// However, for concentric circular geometry, if the circumferantial
// resolution is too coarse, we will end up with some part of the face
// surface that are not taken into account because they do not
// physically overlap any neighbours.  For example, think of 2
// concentric circular patches, slightly rotated one relatively to the
// other.  A good case: the ercoftac conical diffuser, Case0...  GGI
// located right between the cylindrical and conical parts, rotate the
// cylindrical by 15 degrees.  For this case, we will need to devise a
// decent strategy in order to intelligently take care of these
// "missing weights"
//
// The purpose of the ::rescaleWeightingFactors() method is mainly for
// this.
template<class MasterPatch, class SlavePatch>
void newGGIInterpolation<MasterPatch, SlavePatch>::
rescaleWeightingFactors() const
{
    scalarListList& maW = *masterWeightsPtr_;
    scalarListList& saW = *slaveWeightsPtr_;

    // Memorize the largest correction needed in order to provide some
    // basic info to the user
    scalar largestSWC = 0;
    scalar sumSWC = 0;
    scalar curSWC = 0;

    scalar largestMWC = 0;
    scalar sumMWC = 0;
    scalar curMWC = 0;

    // Rescaling the slave weights
    if (debug)
    {
        if
        (
            uncoveredMasterFaces().size() > 0
         || uncoveredSlaveFaces().size() > 0
        )
        {
            InfoIn
            (
                "void newGGIInterpolation<MasterPatch, SlavePatch>::"
                "rescaleWeightingFactors() const"
            )   << "Uncovered faces found.  On master: "
                << uncoveredMasterFaces().size()
                << " on slave: " << uncoveredSlaveFaces().size() << endl;
        }
      }

    forAll (saW, saWi)
    {
        scalar slaveWeightSum = Foam::sum(saW[saWi]);

        if (saW[saWi].size() > 0)
        {
            saW[saWi] = saW[saWi]/slaveWeightSum;

            // Some book-keeping
            curSWC = mag(1.0 - slaveWeightSum);
            largestSWC = Foam::max(largestSWC, curSWC);

            sumSWC += curSWC;
        }
    }

    // Rescaling the master weights
    forAll (maW, maWi)
    {
        scalar masterWeightSum = Foam::sum(maW[maWi]);

        if (maW[maWi].size() > 0)
        {
            maW[maWi] = maW[maWi]/masterWeightSum;

            // Some book-keeping
            curMWC = mag(1.0 - masterWeightSum);
            largestMWC = Foam::max(largestMWC, curMWC);

            sumMWC += curMWC;
        }
    }

    if (debug)
    {
        if (saW.size() > 0 && maW.size() > 0)
        {
            Info<< "  Largest slave weighting factor correction : "
                << largestSWC
                << " average: " << sumSWC/saW.size() << nl
                << "  Largest master weighting factor correction: "
                << largestMWC
                << " average: " << sumMWC/maW.size() << endl;
        }
    }
}


// Find non-overlapping faces from both master and slave patches
// The default non-overlapping criteria is total absence of neighbours.
// Later on, ths criteria will be based on minimum surface intersection, or
// minimum weight factor
template<class MasterPatch, class SlavePatch>
tmp<labelField>
newGGIInterpolation<MasterPatch, SlavePatch>::findNonOverlappingFaces
(
    const scalarListList& patchWeights,
    const scalar& nonOverlapFaceTol   //  = min sum of the neighbour weights
) const
{
    tmp<labelField> tpatchFaceNonOverlapAddr(new labelField());
    labelField& patchFaceNonOverlapAddr = tpatchFaceNonOverlapAddr();

    DynamicList<label, 64> patchFaceNonOverlap(patchWeights.size());

    // Scan the list of patch weights, looking for empty lists
    forAll (patchWeights, paWi)
    {
        scalar sumWeightsFace = sum(patchWeights[paWi]);

        if (sumWeightsFace <= nonOverlapFaceTol)
        {
            // Store local index.
            patchFaceNonOverlap.append(paWi);
        }
    }

    if (patchFaceNonOverlap.size() > 0)
    {
        patchFaceNonOverlapAddr.transfer(patchFaceNonOverlap.shrink());
    }

    if (debug)
    {
        InfoIn("newGGIInterpolation::findNonOverlappingFaces")
            << "   : Found " << patchFaceNonOverlapAddr.size()
            << " non-overlapping faces for this GGI patch" << endl;
    }

    return tpatchFaceNonOverlapAddr;
}


template<class FromPatch, class ToPatch>
void newGGIInterpolation<FromPatch, ToPatch>::
setSourceToTargetPointAddressing
(
    List<labelPair>& sourcePointAddr,   // source point addressing
    scalarField& sourcePointDist,       // source point distances
    vectorField& sourcePointDistVecs,   // source point distance vectors
    const FromPatch& sourcePatch,       // source patch
    const ToPatch& targetPatch,         // target patch
    const labelListList& sourceFaceAddr, // source patch face addressing
    const List< Map<label> >& targetEdgeLoopsMap // target edge loops map
) const
{
    // Procedure
    // for all source points
    //     count how many target faces the source point projects
    //     if the source point (pt S) projects to only one target face
    //         then the projection is unique
    //         the direction is the target face normal direction
    //         the distance is the projected distance to the target face in this
    //         direction
    //     else if the source point (pt S) projects to two target faces
    //         no unique projection
    //         if the two target faces share an edge
    //             find the target edge (edg E) shared by the two target faces
    //             pt S projects to pt S1 on target face 1 and to pt S2 on
    //             target face 2
    //             calculate the perpendicular distance (d1, d2) to edg E from
    //             the projection points (S1, S2)
    //             the direction is the weighted average of the target face
    //             normals, where the weights are based on d1, d2
    //             the distance is taken in this direction; it will project to
    //             one of the faces or the shared edge
    //        else
    //            project to the closest face
    //        end
    //     else if the source point (pt S) projects to greater than two target
    //         faces
    //         no unique projection
    //         if all the target faces share a common point
    //             pt S projects to pt Si on target face i
    //             calculate the distance (di) to pt M from the projection
    //             points (Si)
    //             the direction is the weighted average of the target face
    //             normals, where the weights are based on di
    //             the distance is taken in this direction; it will project to
    //             one of the faces or the shared point
    //         else
    //             project to the closest face
    //         end
    //     else
    //         projection is undefined
    //         if the source point (pt S) projects to any of the target
    //         edges connected to pt M
    //             if the point projects to more than one edge then find the
    //             closest edge
    //             the direction is the normal direction to the edge
    //             the distance is the projected distance to the edge
    //             this direction
    //         else
    //             the direction is from pt S to pt M
    //             the distance is the projected distance in this direction
    //         end
    //     end
    // end

    // Define references for convenience and efficiency
    const labelListList& sourcePointFaces = sourcePatch.pointFaces();
    const faceList& targetFaces = targetPatch.localFaces();
    const pointField& targetPoints = targetPatch.localPoints();
    const pointField& sourcePoints = sourcePatch.localPoints();
    const labelListList& targetPointEdges = targetPatch.pointEdges();
    const edgeList& targetEdges = targetPatch.edges();
    const vectorField targetPointNormals = targetPatch.pointNormals();

    // Allow projection of source points to points/edges on the boundary of the
    // target patch: if true, then all points with possibleMasterFaces will
    // project to the target; if false, then some points may not project to
    // anything
    const labelListList& targetEdgeLoops = targetPatch.edgeLoops();

    forAll(sourcePointAddr, pointI)
    {
        const point& P = sourcePoints[pointI];

        // Find target faces that could possibly be in contact with this point
        // We will use the face addressing for this
        const labelList possibleTargetFaces =
            this->possibleTargetFaces(sourcePointFaces[pointI], sourceFaceAddr);

        if (possibleTargetFaces.size() == 0)
        {
            // We will not calculate the point distances for points that are far
            // from contact i.e. their adjacent faces do not project to faces on
            // the target patch
            continue;
        }

        // Next we will count how many target faces (triangular sub-faces) to
        // which the current point can be projected
        // projectedTriangularSubFaces is a list of triangular sub-faces to
        // which the current point can be projected; each labelPair entry
        // returns the target face index (first) and the face local point index
        // indicating the sub-triangle
        List<labelPair> triSubFacesIndex(0);
        List<vector> triSubFacesProjDirections(0);
        PtrList< triangle<point, point> > triSubFaces(0);
        projectPointToTriangularSubFaces
        (
            triSubFacesIndex,
            triSubFacesProjDirections,
            triSubFaces,
            P,
            possibleTargetFaces,
            targetFaces,
            targetPoints
        );

        if (triSubFacesIndex.size() == 1)
        {
            // Calculate the triangle unit normal
            vector triNormal = triSubFaces[0].normal();
            triNormal /= mag(triNormal);

            // The distance is negative when the point is in contact
            sourcePointDist[pointI] = -triNormal & triSubFacesProjDirections[0];
            sourcePointDistVecs[pointI] = triSubFacesProjDirections[0];

            // The point projects uniquely to one face
            sourcePointAddr[pointI] = triSubFacesIndex[0];
        }
        else if (triSubFacesIndex.size() == 2)
        {
            // Point projects to two faces
            // If they share an edge then we will blend the return
            // directions, else we will pick the face with the smallest
            // distance

            // Triangle A
            const triangle<point, point>& triA = triSubFaces[0];

            // Triangle B
            const triangle<point, point>& triB = triSubFaces[1];

            // Find shared points, if any, of the two triangles
            List<vector> sharedPoints;
            triangleSharedPoints(sharedPoints, triA, triB);

            // If the two sub-faces share an edge
            if (sharedPoints.size() == 2)
            {
                // Find the projection points on both sub-faces
                const point S0 = P + triSubFacesProjDirections[0];
                const point S1 = P + triSubFacesProjDirections[1];

                // Calculate the edge vector direction
                const point& e0 = sharedPoints[0];
                const point& e1 = sharedPoints[1];
                vector e = e1 - e0;
                e /= mag(e);

                // Vectors from S0/S1 to the start of edge e
                const vector S0e0 = e0 - S0;
                const vector S1e0 = e0 - S1;

                // Calculate the perpendicular distance from S0/S1 to edge e
                const scalar d0 = mag((I - sqr(e)) & S0e0);
                const scalar d1 = mag((I - sqr(e)) & S1e0);

                // The projection direction will be a weighted average of
                // the inverted normals of both triangles, where the weights
                // are a function of the perpendicular distances, d0 and d1
                const scalar w0 = d0/(d0 + d1);
                //const scalar w1 = d1/(d0 + d1);
                const scalar w1 = 1.0 - w0;

                // Calculate the triangle unit normals
                vector triN0 = triSubFaces[0].normal();
                triN0 /= mag(triN0);
                vector triN1 = triSubFaces[1].normal();
                triN1 /= mag(triN1);

                // Projection direction is the weighed average of the triangle
                // unit normals
                const vector weightedNormal = w0*triN0 + w1*triN1;
                // We need to be careful with the sign:
                // If eP is in the same direction as weightedNormal then we
                // project in the weightedNormal direction
                // Unit vector from the edge mid-point to point P
                vector eP = P - 0.5*(e0 + e1);
                eP /= mag(eP);
                vector dir = -sign(weightedNormal & eP)*weightedNormal;
                dir /= mag(dir);

                // Calculate the distance along dir to the face
                // We are performing a non-orthogonal projection to the
                // plane of the face
                // e.g. https://computergraphics.stackexchange.com/
                // questions/4360/how-to-project-a-3d-point-onto-a-
                // plane-along-another-axis-vector/4362
                if (w0 > w1)
                {
                    // Point is projected to triangle 0
                    sourcePointDist[pointI] =
                        sign(weightedNormal & eP)
                        *mag(triSubFacesProjDirections[0]/(triN0 & dir));

                    sourcePointAddr[pointI] = triSubFacesIndex[0];
                }
                else
                {
                    // Point is projected to triangle 1
                    sourcePointDist[pointI] =
                        sign(weightedNormal & eP)*
                        mag(triSubFacesProjDirections[1]/(triN1 & dir));

                    sourcePointAddr[pointI] = triSubFacesIndex[1];
                }

                sourcePointDistVecs[pointI] = dir*mag(sourcePointDist[pointI]);
            }
            else
            {
                // The two triangles don't share an edge so we will select
                // the triangle with the shortest projection distance

                // Can we do something better/smoother here?
                // We could check if the triangles share a point and calculate
                // the weighted projection direction based on the distance to
                // this shared point...?

                label closestTriID = 0;
                if
                (
                    mag(triSubFacesProjDirections[1])
                  < mag(triSubFacesProjDirections[0])
                )
                {
                    closestTriID = 1;
                }

                // Calculate triangle unit normal
                vector triNormal = triSubFaces[closestTriID].normal();
                triNormal /= mag(triNormal);

                // Set the distance
                sourcePointDist[pointI] =
                    -triNormal & triSubFacesProjDirections[closestTriID];
                sourcePointDistVecs[pointI] =
                    triSubFacesProjDirections[closestTriID];

                // The point projects to one face
                sourcePointAddr[pointI] = triSubFacesIndex[closestTriID];
            }
        }
        else if (triSubFacesIndex.size() > 2)
        {
            // Point projects to more than two faces

            // Count the number of points common to all the triangles
            List< List<vector> > allSharedPoints(triSubFacesIndex.size() - 1);
            const triangle<point, point>& tri0 = triSubFaces[0];
            for (int triI = 1; triI < triSubFacesIndex.size() - 1; triI++)
            {
                const triangle<point, point>& triangleI = triSubFaces[triI];
                triangleSharedPoints(allSharedPoints[triI], tri0, triangleI);

                if (allSharedPoints[triI].size() == 0)
                {
                    // There are no common shared points
                    break;
                }
            }

            // Check for common shared point
            List<vector> sharedPoints;
            for (int ptI = 0; ptI < allSharedPoints[0].size(); ptI++)
            {
                vector& candidateSharedPoint = allSharedPoints[0][ptI];
                label nTriWithSharedPoint = 1;
                for (int triI = 1; triI < allSharedPoints.size(); triI++)
                {
                    const List<vector>& curSharedPoints = allSharedPoints[triI];
                    bool foundPoint = false;
                    for (int neiI = 0; neiI < curSharedPoints.size(); neiI++)
                    {
                        if (candidateSharedPoint == curSharedPoints[neiI])
                        {
                            foundPoint = true;
                            nTriWithSharedPoint++;
                            break;
                        }
                    }

                    if (!foundPoint)
                    {
                        break;
                    }
                }

                if (nTriWithSharedPoint == allSharedPoints.size())
                {
                    sharedPoints.setSize(1);
                    sharedPoints[0] = candidateSharedPoint;
                }
            }

            // If all the sub-faces share a common point
            if (sharedPoints.size() == 1)
            {
                // If they share a point then we can set the return direction to
                // be the weighted average of the triangle normals, where the
                // weights depend on the distance from the projection point (on
                // each face) to the shared point

                // Shared point
                const vector& M = sharedPoints[0];

                List<point> Si(triSubFacesIndex.size());
                Field<scalar> Wi(triSubFacesIndex.size());
                List<vector> triNi(triSubFacesIndex.size());
                scalar sumW = 0;
                vector weightedNormal = vector::zero;
                forAll(triSubFacesIndex, triI)
                {
                    // Find the projection points on each sub-faces
                    Si[triI] = P + triSubFacesProjDirections[triI];

                    // Calculate the distance from the Si points to point M
                    Wi[triI] = mag(M - Si[triI]);
                    sumW += Wi[triI];

                    // Calculate the triangle unit normals
                    triNi[triI] = triSubFaces[triI].normal();
                    triNi[triI] /= mag(triNi[triI]);

                    // Calculate the weighted normal
                    weightedNormal += Wi[triI]*triNi[triI];
                }

                // Normalise the weights
                Wi /= sumW;
                weightedNormal /= sumW;

                // Unit vector from the opint M to point P
                vector MP = P - M;
                MP /= mag(MP);

                // We need to be careful with the sign:
                // If MP is in the same direction as weightedNormal then we
                // project in the weightedNormal direction
                vector dir = -sign(weightedNormal & MP)*weightedNormal;
                dir /= mag(dir);

                // Find the largest Wi to know which face to project to
                scalar maxWi = 0.0;
                label triWithMaxWi = -1;
                forAll(Wi, triI)
                {
                    if (Wi[triI] > maxWi)
                    {
                        maxWi = Wi[triI];
                        triWithMaxWi = triI;
                    }
                }

                // Point is projected to triangle triWithMaxWi
                sourcePointDist[pointI] =
                    sign(weightedNormal & MP)
                   *mag
                    (
                        triSubFacesProjDirections[triWithMaxWi]
                       /(triNi[triWithMaxWi] & dir)
                    );

                sourcePointAddr[pointI] = triSubFacesIndex[triWithMaxWi];

                sourcePointDistVecs[pointI] = dir*mag(sourcePointDist[pointI]);
            }
            else
            {
                // Project to the closest triangle

                // Like abive, can we do something better/smoother here?

                label closestTriID = -1;
                scalar closestDist = GREAT;

                forAll(triSubFacesIndex, triI)
                {
                    const scalar distI = mag(triSubFacesProjDirections[triI]);

                    if (distI < closestDist)
                    {
                        closestTriID = triI;
                        closestDist = distI;
                    }
                }

                // Calculate triangle unit normal
                vector triNormal = triSubFaces[closestTriID].normal();
                triNormal /= mag(triNormal);

                // Set the distance
                sourcePointDist[pointI] =
                    -triNormal & triSubFacesProjDirections[closestTriID];
                sourcePointDistVecs[pointI] =
                    triSubFacesProjDirections[closestTriID];

                // The point projects to one face
                sourcePointAddr[pointI] = triSubFacesIndex[closestTriID];
            }
        }
        else
        {
            // No point-to-face projection is defined
            // So we will project to an edge or point

            // Find closest vertex on the target
            const label closestTargetPointID =
                findClosestVertex
                (
                    P, possibleTargetFaces, targetFaces, targetPoints
                );
            const point& M =
                targetPoints[closestTargetPointID];

            // For all edges connected to the target point, we will check if
            // the source point projects to the edge

            // Get the edges attached to the target point
            const labelList& closestEdges =
                targetPointEdges[closestTargetPointID];

            // Find the closest edge to project to, if any
            label closestEdgeI = -1;
            scalar minDist = GREAT;
            forAll(closestEdges, eI)
            {
                const label curEdgeID = closestEdges[eI];
                const edge& curEdge = targetEdges[curEdgeID];

                // Vector from M to P
                const vector MP = P - M;

                // Vertex at the other end of the edge
                const label otherPointID =
                    curEdge.otherVertex(closestTargetPointID);

                // Unit edge vector away from M
                vector e = targetPoints[otherPointID] - M;
                const scalar magE = mag(e);
                e /= magE;

                // Project point P to the edge and calculate the distance
                // along the edge to the projection point
                const scalar d = e & MP;

                // If 0 <= d <= magE, then P projects to the edge
                if (d >= 0 && d <= magE)
                {
                    const scalar dist = mag((I - sqr(e)) & MP);

                    if (dist < minDist)
                    {
                        minDist = dist;
                        closestEdgeI = eI;
                    }
                }
            }

            // If the point projects to an edge then perform the projection,
            // else we will project to the point M
            if (closestEdgeI > -1)
            {
                const label curEdgeID = closestEdges[closestEdgeI];
                const edge& curEdge = targetEdges[curEdgeID];

                // Vector from M to P
                const vector MP = P - M;

                // Vertex at the other end of the edge
                const label otherPointID =
                    curEdge.otherVertex(closestTargetPointID);

                // Unit edge vector away from M
                vector e = targetPoints[otherPointID] - M;
                const scalar magE = mag(e);
                e /= magE;

                // To check if the distance is positive or negative, we
                // will compare the projection vector with the edge normal
                // vector, where we will define the edge normal vector as the
                // average of the edge-start and edge-end point normals
                const vector& startNormal = targetPointNormals[curEdge.start()];
                const vector& endNormal = targetPointNormals[curEdge.end()];
                vector edgeN = 0.5*(startNormal + endNormal);
                const scalar magEdgeN = mag(edgeN);
                if (magEdgeN < SMALL)
                {
                    FatalErrorIn
                    (
                        "template<class FromPatch, class ToPatch>\n"
                        "void newGGIInterpolation<FromPatch, ToPatch>::\n"
                        "setSourceToTargetPointAddressing\n"
                        "(\n"
                        "    List<labelPair>& sourcePointAddressing,\n"
                        "    scalarField& sourcePointDistance,\n"
                        "    vectorField& sourcePointDistVecs,\n"
                        "    const FromPatch& sourcePatch,\n"
                        "    const ToPatch& targetPatch,\n"
                        "    const labelListList& sourceAddr\n"
                        ") const"
                    )   << "Something went wrong projecting to an edge!"
                        << " The edge start point normal is the neative of the "
                        << "end point normal!"
                        << abort(FatalError);
                }
                edgeN /= mag(edgeN);

                // Projection vector from point P to the closest point
                // on the edge
                sourcePointDist[pointI] =
                    sign(edgeN & MP)*mag((I - sqr(e)) & MP);
                sourcePointDistVecs[pointI] = -(I - sqr(e)) & MP;

                // The point projects to one an edge so there is no
                // triangle address to store here, so we will leave it
                // as -1
                sourcePointAddr[pointI] = labelPair(-1, -1);

                if (!projectPointsToPatchBoundary_)
                {
                    correctBoundaryPointProjection
                    (
                        targetPatch,
                        closestTargetPointID,
                        targetEdgeLoops,
                        targetEdgeLoopsMap,
                        MP,
                        sourcePointDist[pointI],
                        sourcePointDistVecs[pointI]
                    );
                }
            }
            else
            {
                // If the point does not project to any edge or projects to more
                // than one edge then we will project directly to the point M

                // Vector from M to P
                const vector MP = P - M;

                // Get the point normal
                const vector& pointN = targetPointNormals[closestTargetPointID];

                // The distance is negative when the point is in contact
                sourcePointDist[pointI] = sign(pointN & MP)*mag(MP);
                sourcePointDistVecs[pointI] = -MP;

                // The point projects to a point so there is no triangle
                // address to store here; we will leave it as -1
                sourcePointAddr[pointI] = labelPair(-1, -1);

                if (!projectPointsToPatchBoundary_)
                {
                    correctBoundaryPointProjection
                    (
                        targetPatch,
                        closestTargetPointID,
                        targetEdgeLoops,
                        targetEdgeLoopsMap,
                        MP,
                        sourcePointDist[pointI],
                        sourcePointDistVecs[pointI]
                    );
                }
            }
        }
    }

    // Check orientation
    if (checkPointDistanceOrientations_)
    {
        const pointField& sourcePointNormals = sourcePatch.pointNormals();
        const vectorField& targetFaceNormals = targetPatch.faceNormals();

        scalarField orientation(sourcePointAddr.size(), 0);

        label nIncorrectPoints = 0;

        forAll(sourcePointAddr, pointI)
        {
            if (sourcePointAddr[pointI].first() != -1)
            {
                orientation[pointI] =
                    (
                        sourcePointNormals[pointI]
                      & targetFaceNormals[sourcePointAddr[pointI].first()]
                    );

                if (orientation[pointI] > -SMALL)
                {
                    nIncorrectPoints++;

                    // Ignore these points
                    sourcePointAddr[pointI] = labelPair(-1, -1);
                    sourcePointDist[pointI] = GREAT;
                }
            }
        }
    }

    if (debug)
    {
        Info<< "GGI, source point distance"
            << ", max: " << max(sourcePointDist)
            << ", min: " << min(sourcePointDist)
            << ", avg: " << average(sourcePointDist) << endl;
    }
}


template<class FromPatch, class ToPatch>
void newGGIInterpolation<FromPatch, ToPatch>::
setSourceToTargetPointAddressingPrev
(
    List<labelPair>& sourcePointAddr,   // source point addressing
    scalarField& sourcePointDist,       // source point distances
    const FromPatch& sourcePatch,       // source patch
    const ToPatch& targetPatch,         // target patch
    const labelListList& sourceFaceAddr // source patch face addressing
) const
{
    // Take references for convenience and efficiency
    const labelListList& pointFaces = sourcePatch.pointFaces();
    const faceList& targetFaces = targetPatch.localFaces();
    const pointField& targetPoints = targetPatch.localPoints();
    const pointField& sourcePoints = sourcePatch.localPoints();

    forAll(sourcePointAddr, pointI)
    {
        const point& P = sourcePoints[pointI];

        labelHashSet possibleTargetFacesSet;

        const labelList& curPointFaces = pointFaces[pointI];
        forAll(curPointFaces, faceI)
        {
            const label curFace = curPointFaces[faceI];

            const labelList& curTargetFaces = sourceFaceAddr[curFace];

            forAll(curTargetFaces, fI)
            {
                if (!possibleTargetFacesSet.found(curTargetFaces[fI]))
                {
                    possibleTargetFacesSet.insert(curTargetFaces[fI]);
                }
            }
        }

        labelList possibleTargetFaces = possibleTargetFacesSet.toc();

        scalar MinEta = -GREAT;
        labelPair faceTriangle(-1, -1);
        scalar distance = GREAT;

        forAll(possibleTargetFaces, faceI)
        {
            const label curTargetFace = possibleTargetFaces[faceI];

            const face& f = targetFaces[curTargetFace];

            const point ctr = Foam::average(f.points(targetPoints));

            //point nextPoint = ctr;

            for (label pI = 0; pI < f.size(); pI++)
            {
                const point nextPoint = targetPoints[f.nextLabel(pI)];

                triPointRef t
                (
                    targetPoints[f[pI]],
                    nextPoint,
                    ctr
                );

                vector n = t.normal();
                const scalar A = mag(n);
                n /= A;

                // Intersection point
                const point I = P + n*(n & (t.a() - P));

                // Areal coordinates
                scalarField eta(3, 0);

                eta[0] = (triPointRef(I, t.b(), t.c()).normal() & n)/A;
                eta[1] = (triPointRef(I, t.c(), t.a()).normal() & n)/A;
                eta[2] = (triPointRef(I, t.a(), t.b()).normal() & n)/A;

                const scalar minEta = min(eta);

                // PC/PDJ/TT: 28-11-17: for projections within the face, the
                // sum of mag(eta[0]) + mag(eta[1]) + mag(eta[2]) must be equal
                // to one if the point is projected within the triangle
                // We will allow 5% to account for convex corners, but we must
                // come back to solve these convex corners correctly!!!
                // PC: 01-11-18: I'm back: the default behaviour is more robust
                // so we will not force the point to be projected within the
                // face
                // if (mag(eta[0]) + mag(eta[1]) + mag(eta[2]) < 1.05)
                {
                    if (minEta > MinEta)
                    {
                        MinEta = minEta;
                        faceTriangle.first() = curTargetFace;
                        faceTriangle.second() = pI;

                        distance = ((P - I) & n);
                    }
                }
            }
        }

        sourcePointAddr[pointI] = faceTriangle;
        sourcePointDist[pointI] = distance;
    }

    if (debug)
    {
        Info<< "GGI, source point distance"
            << ", max: " << max(sourcePointDist)
            << ", min: " << min(sourcePointDist)
            << ", avg: " << average(sourcePointDist) << endl;
    }

    // Check orientation
    if (checkPointDistanceOrientations_)
    {
        const pointField& sourcePointNormals =
            sourcePatch.pointNormals();

        const vectorField& targetFaceNormals =
            targetPatch.faceNormals();

        scalarField orientation(sourcePointAddr.size(), 0);

        label nIncorrectPoints = 0;

        forAll(sourcePointAddr, pointI)
        {
            if (sourcePointAddr[pointI].first() != -1)
            {
                orientation[pointI] =
                    (
                        sourcePointNormals[pointI]
                      & targetFaceNormals[sourcePointAddr[pointI].first()]
                    );

                if (orientation[pointI] > -SMALL)
                {
                    nIncorrectPoints++;

                    // Ignore these points
                    sourcePointAddr[pointI] = labelPair(-1, -1);
                    sourcePointDist[pointI] = GREAT;
                }
            }
        }
    }
}


template<class FromPatch, class ToPatch>
tmp<labelField> newGGIInterpolation<FromPatch, ToPatch>::possibleTargetFaces
(
    const labelList& sourceFaces,
    const labelListList& sourceFaceAddr
) const
{
    tmp<labelField> tpossibleFaces(new labelField());
    labelField& possibleFaces = tpossibleFaces();

    labelHashSet possibleFacesSet;

    forAll(sourceFaces, faceI)
    {
        const label curFace = sourceFaces[faceI];
        const labelList& curTargetFaces = sourceFaceAddr[curFace];

        forAll(curTargetFaces, fI)
        {
            if (!possibleFacesSet.found(curTargetFaces[fI]))
            {
                possibleFacesSet.insert(curTargetFaces[fI]);
            }
        }
    }

    possibleFaces = possibleFacesSet.toc();

    return tpossibleFaces;
}

template<class FromPatch, class ToPatch>
void newGGIInterpolation<FromPatch, ToPatch>::
projectPointToTriangularSubFaces
(
    List<labelPair>& triSubFacesIndex,
    List<vector>& triSubFacesProjDirections,
    PtrList< triangle<point, point> >& triSubFaces,
    const point& P,
    const labelList& possibleTargetFaces,
    const faceList& targetFaces,
    const pointField& targetPoints
) const
{
    DynamicList<labelPair> faceTrianglesIndex(10);
    DynamicList<vector> faceTrianglesDistVecs(10);

    forAll(possibleTargetFaces, faceI)
    {
        const label curFaceID = possibleTargetFaces[faceI];
        const face& curFace = targetFaces[curFaceID];

        // Centre of the face: average position of the points
        const point ctr = Foam::average(curFace.points(targetPoints));

        //point nextPoint = ctr;

        for (label pI = 0; pI < curFace.size(); pI++)
        {
            // Next point on the face following the right-hand rule
            const point nextPoint = targetPoints[curFace.nextLabel(pI)];

            // Create the triangular sub-face
            const triPointRef t
            (
                targetPoints[curFace[pI]],
                nextPoint,
                ctr
            );

            // Area vector of the triangle (we make this a unit normal below)
            vector n = t.normal();

            // Area magnitude of the triangle
            const scalar A = mag(n);

            // Unit normal of the triangle
            n /= A;

            // Intersection point when P is projected in the normal direction to
            // the triangle
            const point I = P + n*(n & (t.a() - P));

            // Areal coordinates within the triangle e.g. see
            // https://en.wikipedia.org/wiki/Barycentric_coordinate_system
            scalarField eta(3, 0);
            eta[0] = (triPointRef(I, t.b(), t.c()).normal() & n)/A;
            eta[1] = (triPointRef(I, t.c(), t.a()).normal() & n)/A;
            eta[2] = (triPointRef(I, t.a(), t.b()).normal() & n)/A;

            // If I is within the face, the min(eta) must be greater than or
            // equal to 0 and the sum (mag(eta[0]) + mag(eta[1]) + mag(eta[2]))
            // must be equal to 1.
            if (min(eta) > 0.0)
            {
                faceTrianglesIndex.append(labelPair(curFaceID, pI));

                // The distance vector from point P to the point I on the face
                faceTrianglesDistVecs.append(sqr(n) & (I - P));

                // The distance vector is given as
                const label oldSize = triSubFaces.size();
                triSubFaces.resize(oldSize + 1);
                triSubFaces.set
                (
                    oldSize,
                    new triangle<point, point>(t.a(), t.b(), t.c())
                );
            }
        }
    }

    triSubFacesIndex = faceTrianglesIndex;
    triSubFacesProjDirections = faceTrianglesDistVecs;
}


template<class FromPatch, class ToPatch>
void newGGIInterpolation<FromPatch, ToPatch>::
triangleSharedPoints
(
    List<point>& sharedPoints,
    const triangle<point, point>& triA,
    const triangle<point, point>& triB
) const
{
    int nSharedPoints = 0;
    List<vector> sp(3, vector::zero);

    if
    (
        triA.a() == triB.a() || triA.a() == triB.b() || triA.a() == triB.c()
    )
    {
        sp[nSharedPoints] = triA.a();
        nSharedPoints++;
    }

    if
    (
        triA.b() == triB.a() || triA.b() == triB.b() || triA.b() == triB.c()
    )
    {
        sp[nSharedPoints] = triA.b();
        nSharedPoints++;
    }

    if
    (
        triA.c() == triB.a() || triA.c() == triB.b() || triA.c() == triB.c()
    )
    {
        sp[nSharedPoints] = triA.c();
        nSharedPoints++;
    }

    sharedPoints.resize(nSharedPoints);
    forAll(sharedPoints, pI)
    {
        sharedPoints[pI] = sp[pI];
    }
}


template<class MasterPatch, class SlavePatch>
void newGGIInterpolation<MasterPatch, SlavePatch>::
correctBoundaryPointProjection
(
    const MasterPatch& targetPatch,
    const label closestTargetPointID,
    const labelListList& targetEdgeLoops,
    const List< Map<label> >& targetEdgeLoopsMap,
    const vector& MP,
    scalar& sourcePointDist,
    vector& sourcePointDistVecs
) const
{
    forAll(targetEdgeLoopsMap, loopI)
    {
        const Map<label>& targetEdgeLoopsMapI =
            targetEdgeLoopsMap[loopI];

        if (targetEdgeLoopsMapI.found(closestTargetPointID))
        {
            // M is on the boundary of the target patch

            // Calculate the vector pointing out from the patch
            // interior, that is perpendicular to the point
            // normal and to the average direction vectors of
            // the adjacent boundary edges. This is equivalent
            // to the "m" vector in the finite area method
            vector m = vector::zero;
            boundaryPointBiNormal
            (
                m,
                targetPatch,
                closestTargetPointID,
                targetEdgeLoops[loopI],
                targetEdgeLoopsMapI
            );

            if ((MP & m) > SMALL)
            {
                sourcePointDist = GREAT;
                sourcePointDistVecs = vector::zero;

                return;
            }
        }
    }
}

template<class MasterPatch, class SlavePatch>
void newGGIInterpolation<MasterPatch, SlavePatch>::boundaryPointBiNormal
(
    vector& m,
    const MasterPatch& targetPatch,
    const label pointID,
    const labelList& targetBoundaryPoints,
    const Map<label>& targetBoundaryPointsMap
) const
{
    // Take references
    const labelListList& pointFaces = targetPatch.pointFaces();
    const faceList& faces = targetPatch.localFaces();
    const pointField& points = targetPatch.localPoints();

    // Find index in boundaryPoints list, which is an
    // ordered list of the boundary points
    Map<label>::const_iterator iter = targetBoundaryPointsMap.find(pointID);

    if (iter == Map<label>::end())
    {
        FatalErrorIn
        (
            "template<class MasterPatch, class SlavePatch>\n"
            "vector newGGIInterpolation<MasterPatch, SlavePatch>::\n"
            "boundaryPointBiNormal\n"
            "(\n"
            "    vector& m,\n"
            "    const MasterPatch& targetPatch,\n"
            "    const label pointID,\n"
            "    const labelList& targetBoundaryPoints,\n"
            "    const Map<label> targetBoundaryPointsMap\n"
            ") const"
        )   << "Could not point in the boundary points map!"
            << abort(FatalError);
    }

    const label bpI = iter();

    // Find the next and previous points along the boundary
    label nextBpI = bpI + 1;
    label prevBpI = bpI - 1;
    if (prevBpI == -1)
    {
        prevBpI = targetBoundaryPoints.size() - 1;
    }
    if (nextBpI == targetBoundaryPoints.size())
    {
        nextBpI = 0;
    }

    const label nextPointID = targetBoundaryPoints[nextBpI];
    const label prevPointID = targetBoundaryPoints[prevBpI];

    // Find the face containing bpI and prevBpI
    const labelList& curPointFaces = pointFaces[pointID];
    label nextFaceID = -1;
    label prevFaceID = -1;
    forAll(curPointFaces, pI)
    {
        const label faceID = curPointFaces[pI];
        const face& curFace = faces[faceID];

        forAll(curFace, fpI)
        {
            if (curFace[fpI] == prevPointID)
            {
                prevFaceID = faceID;
                break;
            }
            else if (curFace[fpI] == nextPointID)
            {
                nextFaceID = faceID;
                break;
            }
        }
    }

    if (prevFaceID == -1 || nextFaceID == -1)
    {
        FatalErrorIn
        (
            "template<class MasterPatch, class SlavePatch>\n"
            "vector newGGIInterpolation<MasterPatch, SlavePatch>::\n"
            "boundaryPointBiNormal\n"
            "(\n"
            "    vector& m,\n"
            "    const MasterPatch& targetPatch,\n"
            "    const label pointID,\n"
            "    const labelList& targetBoundaryPoints,\n"
            "    const Map<label> targetBoundaryPointsMap\n"
            ") const"
        )   << "Could not find previous and/or next faces at a boundary point!"
            << nl << "PrevPoint: " << points[prevPointID]
            << ", Point: " << points[bpI] << ", NextPoint: "
            << points[nextPointID] << nl
            << "prevBpI: " << prevBpI << ", bpI: " << bpI
            << ", nextBpI: " << nextBpI
            << nl << "prevFaceID: " << prevFaceID
            << ", nextFaceID: " << nextFaceID
            << abort(FatalError);
    }

    const face& prevFace = faces[prevFaceID];
    const face& nextFace = faces[nextFaceID];

    // Find vector from prevFace/nextFace centre to point M
    const point& M = points[pointID];
    const vector prevM = M - prevFace.centre(points);
    const vector nextM = M - nextFace.centre(points);

    // We will take the average of prevM and nextM as the biNormal direction
    m = 0.5*(prevM + nextM);
    m /= mag(m);
}


template<class FromPatch, class ToPatch>
label newGGIInterpolation<FromPatch, ToPatch>::
findClosestVertex
(
    const point& P,
    const labelList& possibleTargetFaces,
    const faceList& targetFaces,
    const pointField& targetPoints
) const
{
    scalar minDist = GREAT;
    label closestPointID = -1;

    forAll(possibleTargetFaces, faceI)
    {
        const label curFaceID = possibleTargetFaces[faceI];
        const face& curFace = targetFaces[curFaceID];

        forAll(curFace, pI)
        {
            const label curPointID = curFace[pI];
            const point& curPoint = targetPoints[curPointID];
            const scalar dist = mag(P - curPoint);

            if (dist < minDist)
            {
                minDist = dist;
                closestPointID = curPointID;
            }
        }
    }

    if (closestPointID == -1)
    {
        FatalErrorIn("newGGIInterpolation::findClosestVertex(...)")
            << "No point was found!" << abort(FatalError);
    }

    return closestPointID;
}


template<class FromPatch, class ToPatch>
void newGGIInterpolation<FromPatch, ToPatch>::
calcMasterPointAddressing() const
{
    if (masterPointAddressingPtr_)
    {
        FatalErrorIn
        (
            "void newGGIInterpolation::calcMasterPointAddressing() const"
        )   << "Pointer already set!" << abort(FatalError);
    }

    masterPointAddressingPtr_ =
        new List<labelPair>
        (
            this->masterPatch().nPoints(),
            labelPair(-1,-1)
        );
    List<labelPair>& masterPointAddr = *masterPointAddressingPtr_;

    masterPointDistancePtr_ =
        new scalarField
        (
            this->masterPatch().nPoints(),
            GREAT
        );
    scalarField& masterPointDist = *masterPointDistancePtr_;

    masterPointDistanceVectorsPtr_ =
        new vectorField
        (
            this->masterPatch().nPoints(),
            vector::zero
        );
    vectorField& masterPointDistVecs = *masterPointDistanceVectorsPtr_;

    if (useNewPointDistanceMethod_)
    {
        // Set the point addressing from the source patch points to the target
        // patch faces
        setSourceToTargetPointAddressing
        (
            masterPointAddr,     // source point addressing
            masterPointDist,     // source point distances
            masterPointDistVecs, // source point distance vectors
            masterPatch(),       // source patch
            slavePatch(),        // target patch
            masterAddr(),        // source patch face addressing
            slaveEdgeLoopsMap()  // target edge loops map
        );
    }
    else
    {
        // Previous method
        // Note: the distance vectors are not set and left as zero
        setSourceToTargetPointAddressingPrev
        (
            masterPointAddr,     // source point addressing
            masterPointDist,     // source point distances
            masterPatch(),       // source patch
            slavePatch(),        // target patch
            masterAddr()         // source patch face addressing
        );
    }
}


template<class FromPatch, class ToPatch>
void newGGIInterpolation<FromPatch, ToPatch>::
calcMasterPointWeights() const
{
    // Find master point weights
    if (masterPointWeightsPtr_)
    {
        FatalErrorIn
        (
            "void newGGIInterpolation::"
            "calcMasterPointAddressing() const"
        )   << "Master point weights already exist"
            << abort(FatalError);
    }

    masterPointWeightsPtr_ =
        new FieldField<Field, scalar>(this->masterPatch().nPoints());
    FieldField<Field, scalar>& masterPointWeights = *masterPointWeightsPtr_;

    const faceList& slaveFaces = this->slavePatch().localFaces();
    const pointField& slavePoints = this->slavePatch().localPoints();

    const pointField& masterPoints = this->masterPatch().localPoints();

    const List<labelPair>& addr = this->masterPointAddr();

    forAll(masterPointWeights, pointI)
    {
        if (addr[pointI].first() != -1)
        {
            const point& P = masterPoints[pointI];

            const face& hitFace =
                slaveFaces[addr[pointI].first()];

            point ctr = Foam::average(hitFace.points(slavePoints));

            label pI = addr[pointI].second();

            triPointRef t
            (
                slavePoints[hitFace[pI]],
                slavePoints[hitFace.nextLabel(pI)],
                ctr
            );

            vector n = t.normal();
            n /= mag(n);

            // Intersection point
            point I = P + n*(n&(t.a() - P));

            masterPointWeights.set(pointI, scalarField(3));

            masterPointWeights[pointI][0] = t.Ni(0, I);
            masterPointWeights[pointI][1] = t.Ni(1, I);
            masterPointWeights[pointI][2] = t.Ni(2, I);
        }
        else
        {
            masterPointWeights.set(pointI, scalarField(0));
        }
    }
}

template<class FromPatch, class ToPatch>
void newGGIInterpolation<FromPatch, ToPatch>::calcSlavePointAddressing() const
{
    if (slavePointAddressingPtr_)
    {
        FatalErrorIn
        (
            "void newGGIInterpolation::calcSlavePointAddressing() const"
        )   << "Pointer already set!" << abort(FatalError);
    }

    slavePointAddressingPtr_ =
        new List<labelPair>
        (
            this->slavePatch().nPoints(),
            labelPair(-1,-1)
        );
    List<labelPair>& slavePointAddr = *slavePointAddressingPtr_;

    slavePointDistancePtr_ =
        new scalarField
        (
            this->slavePatch().nPoints(),
            GREAT
        );
    scalarField& slavePointDist = *slavePointDistancePtr_;

    slavePointDistanceVectorsPtr_ =
        new vectorField
        (
            this->slavePatch().nPoints(),
            vector::zero
        );
    vectorField& slavePointDistVecs = *slavePointDistanceVectorsPtr_;

    if (useNewPointDistanceMethod_)
    {
        // Set the point addressing from the source patch points to the target
        // patch faces
        setSourceToTargetPointAddressing
        (
            slavePointAddr,      // source point addressing
            slavePointDist,      // source point distances
            slavePointDistVecs,  // source point distance vectors
            slavePatch(),        // source patch
            masterPatch(),       // target patch
            slaveAddr(),         // source patch face addressing
            masterEdgeLoopsMap() // target edge loops map
        );
    }
    else
    {
        // Set the point addressing from the source patch points to the target
        // patch faces using the previous method
        // Note: the distance vectors are not set and are left as zero
        setSourceToTargetPointAddressingPrev
        (
            slavePointAddr,     // source point addressing
            slavePointDist,     // source point distances
            slavePatch(),       // source patch
            masterPatch(),      // target patch
            slaveAddr()         // source patch face addressing
        );
    }
}

template<class FromPatch, class ToPatch>
void newGGIInterpolation<FromPatch, ToPatch>::
calcSlavePointWeights() const
{
    // Find master point weights
    if (slavePointWeightsPtr_)
    {
        FatalErrorIn
        (
            "void newGGIInterpolation::"
            "calcSlavePointAddressing() const"
        )   << "Slave point weights already exist"
            << abort(FatalError);
    }

    slavePointWeightsPtr_ =
        new FieldField<Field, scalar>(this->slavePatch().nPoints());
    FieldField<Field, scalar>& slavePointWeights = *slavePointWeightsPtr_;

    const faceList& masterFaces = this->masterPatch().localFaces();
    const pointField& masterPoints = this->masterPatch().localPoints();

    const pointField& slavePoints = this->slavePatch().localPoints();

    const List<labelPair>& addr = this->slavePointAddr();

    forAll(slavePointWeights, pointI)
    {
        if (addr[pointI].first() != -1)
        {
            const point& P = slavePoints[pointI];

            const face& hitFace =
                masterFaces[addr[pointI].first()];

            point ctr = Foam::average(hitFace.points(masterPoints));

            label pI = addr[pointI].second();

            triPointRef t
            (
                masterPoints[hitFace[pI]],
                masterPoints[hitFace.nextLabel(pI)],
                ctr
            );

            vector n = t.normal();
            n /= mag(n);

            // Intersection point
            point I = P + n*(n&(t.a() - P));

            slavePointWeights.set(pointI, scalarField(3));

            slavePointWeights[pointI][0] = t.Ni(0, I);
            slavePointWeights[pointI][1] = t.Ni(1, I);
            slavePointWeights[pointI][2] = t.Ni(2, I);
        }
        else
        {
            slavePointWeights.set(pointI, scalarField(0));
        }
    }
}


template<class FromPatch, class ToPatch>
void newGGIInterpolation<FromPatch, ToPatch>::calcMasterEdgeLoopsMap() const
{
    if (masterEdgeLoopsMap_.size() > 0)
    {
        FatalErrorIn
        (
            "template<class FromPatch, class ToPatch>"
            "void newGGIInterpolation<FromPatch, ToPatch>::"
            "calcMasterEdgeLoopsMap() const"
        )   << "Pointer already set" << abort(FatalError);
    }

    const labelListList& masterEdgeLoops = masterPatch().edgeLoops();

    masterEdgeLoopsMap_.setSize(masterEdgeLoops.size());

    forAll(masterEdgeLoopsMap_, loopI)
    {
        const labelList& masterEdgeLoopsI = masterEdgeLoops[loopI];
        Map<label>& masterEdgeLoopsMapI = masterEdgeLoopsMap_[loopI];

        forAll(masterEdgeLoopsI, pI)
        {
            masterEdgeLoopsMapI.insert(masterEdgeLoopsI[pI], pI);
        }
    }
}


template<class FromPatch, class ToPatch>
void newGGIInterpolation<FromPatch, ToPatch>::calcSlaveEdgeLoopsMap() const
{
    if (slaveEdgeLoopsMap_.size() > 0)
    {
        FatalErrorIn
        (
            "template<class FromPatch, class ToPatch>"
            "void newGGIInterpolation<FromPatch, ToPatch>::"
            "calcSlaveEdgeLoopsMap() const"
        )   << "Pointer already set" << abort(FatalError);
    }

    const labelListList& slaveEdgeLoops = slavePatch().edgeLoops();

    slaveEdgeLoopsMap_.setSize(slaveEdgeLoops.size());

    forAll(slaveEdgeLoopsMap_, loopI)
    {
        const labelList& slaveEdgeLoopsI = slaveEdgeLoops[loopI];

        forAll(slaveEdgeLoopsI, pI)
        {
            slaveEdgeLoopsMap_[loopI].insert(slaveEdgeLoopsI[pI], pI);
        }
    }
}


template<class FromPatch, class ToPatch>
const List< Map<label> >&
newGGIInterpolation<FromPatch, ToPatch>::masterEdgeLoopsMap() const
{
    if (masterEdgeLoopsMap_.size() == 0)
    {
        calcMasterEdgeLoopsMap();
    }

    return masterEdgeLoopsMap_;
}


template<class FromPatch, class ToPatch>
const List< Map<label> >&
newGGIInterpolation<FromPatch, ToPatch>::slaveEdgeLoopsMap() const
{
    if (slaveEdgeLoopsMap_.size() == 0)
    {
        calcSlaveEdgeLoopsMap();
    }

    return slaveEdgeLoopsMap_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
