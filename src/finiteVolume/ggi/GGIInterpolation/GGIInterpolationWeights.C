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

Description
    Mass-conservative face interpolation of face data between two
    primitivePatches

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

Modification by:
    Martin Beaudoin, Hydro-Quebec, (2008)

\*---------------------------------------------------------------------------*/

#include "GGIInterpolationTemplate.H"
#include "objectHit.H"
#include "boolList.H"
#include "DynamicList.H"
#include "dimensionedConstants.H"
#include "triPointShapeRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// template<class MasterPatch, class SlavePatch>
// const Foam::debug::tolerancesSwitch
// GGIInterpolation<MasterPatch, SlavePatch>::areaErrorTol_
// (
//     "GGIAreaErrorTol",
//     1.0e-8,
//     "Minimum GGI face to face intersection area.  "
//     "The smallest accepted GGI weighting factor."
// );
//
//
// template<class MasterPatch, class SlavePatch>
// const Foam::debug::tolerancesSwitch
// GGIInterpolation<MasterPatch, SlavePatch>::featureCosTol_
// (
//     "GGIFeatureCosTol",
//     0.8,
//     "Minimum cosine value between 2 GGI patch neighbouring facet normals."
// );
//
//
// template<class MasterPatch, class SlavePatch>
// Foam::debug::tolerancesSwitch
// GGIInterpolation<MasterPatch, SlavePatch>::uncoveredFaceAreaTol_
// (
//     "GGIUncoveredFaceAreaTol",
//     0.999,
//     "Fraction of face area mismatch (sum of weights) to consider a face "
//     "as uncovered, i.e. not to rescale weights."
// );
template<class MasterPatch, class SlavePatch>
const Foam::scalar
GGIInterpolation<MasterPatch, SlavePatch>::areaErrorTol_ = 1.0e-8;

template<class MasterPatch, class SlavePatch>
const Foam::scalar
GGIInterpolation<MasterPatch, SlavePatch>::featureCosTol_ = 0.8;


template<class MasterPatch, class SlavePatch>
Foam::scalar
GGIInterpolation<MasterPatch, SlavePatch>::uncoveredFaceAreaTol_ = 0.999;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
void GGIInterpolation<MasterPatch, SlavePatch>::calcAddressing() const
{
    if
    (
        masterAddrPtr_
     || masterWeightsPtr_
     || slaveAddrPtr_
     || slaveWeightsPtr_
     || uncoveredMasterAddrPtr_
     || partiallyUncoveredMasterAddrPtr_
     || masterFaceUncoveredFractionsPtr_
     || uncoveredSlaveAddrPtr_
     || partiallyUncoveredSlaveAddrPtr_
     || slaveFaceUncoveredFractionsPtr_
    )
    {
        FatalErrorIn
        (
            "void GGIInterpolation<MasterPatch, SlavePatch>::"
            "calcAddressing() const"
        )   << "Addressing already calculated"
            << abort(FatalError);
    }

    if (debug)
    {
        InfoIn
        (
            "void GGIInterpolation<MasterPatch, SlavePatch>::"
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

    if (reject_ == AABB)
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
            "void GGIInterpolation<MasterPatch, SlavePatch>::"
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

    const pointField& masterPatchPoints = masterPatch_.points();
    const vectorField masterPatchNormals = masterPatch_.faceNormals();

    // // ZT, 05/07/2014
    // const vectorField& slavePatchNormals = slavePatch_.faceNormals();

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
                "void GGIInterpolation<MasterPatch, SlavePatch>::"
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
                slavePatch_[curCMN[neighbI]].points(slavePatch_.points());

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

            // // ZT, 05/07/2014
            // scalar orientation =
            //     (
            //         masterPatchNormals[faceMi]
            //       & slavePatchNormals[curCMN[neighbI]]
            //     );

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
                    sqrt(areaErrorTol_) // distErrorTol
                )
             // && (orientation < -SMALL) // ZT, 05/07/2014
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

                scalar intersectionTestArea =
                    Foam::max
                    (
                        VSMALL,
                        areaErrorTol_*
                        Foam::max
                        (
                            surfaceAreaMasterPointsInUV,
                            surfaceAreaNeighbPointsInUV
                        )
                    );

                // Fix: previously checked for VSMALL.
                // HJ, 19/Sep/2016
                if (intersectionArea > intersectionTestArea)
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
//                     WarningIn
//                     (
//                         "GGIInterpolation<MasterPatch, SlavePatch>::"
//                         "calcAddressing()"
//                     )   << "polygonIntersection is returning a "
//                         << "zero surface area between " << nl
//                         << "     Master face: " << faceMi
//                         << " and Neighbour face: " << curCMN[neighbI]
//                         << " intersection area = " << intersectionArea << nl
//                         << "Please check the two quick-check algorithms for "
//                         << "GGIInterpolation.  Something is  missing." << endl;
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

        // Current master face index
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

    // Calculate master and slave partially covered addressing

    // Note: This function allocates:
    // 1. partiallyUncoveredMasterAddrPtr_
    // 2. masterFaceUncoveredFractionsPtr_
    calcPartiallyCoveredFaces
    (
        maW,
        masterNonOverlapFaceTol_,
        true // This is master
    );

    // Note: this function allocates:
    // 1. partiallyUncoveredSlaveAddrPtr_
    // 2. slaveFaceUncoveredFractionsPtr_
    calcPartiallyCoveredFaces
    (
        saW,
        slaveNonOverlapFaceTol_,
        false // This is not master
    );


    // Rescaling the weighting factors so they will sum up to 1.0
    // See the comment for the method ::rescaleWeightingFactors() for
    // more information.  By default, we always rescale.  But for some
    // special kind of GGI interpolation, like the mixingPlaneGGI,
    // then we need the brute values, so no rescaling in that
    // case. Hence the little flag rescaleGGIWeightingFactors_

    // Not parallelised.  HJ, 27/Apr/2016. Rescaling is not performed for
    // partially overlapping faces for their correct treatment. VV, 16/Oct/2017.
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
void GGIInterpolation<MasterPatch, SlavePatch>::rescaleWeightingFactors() const
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

    // Note: do not rescale weighting factors for partially covered faces
    if (!partiallyUncoveredMasterAddrPtr_ || !partiallyUncoveredSlaveAddrPtr_)
    {
        FatalErrorIn
        (
            "void GGIInterpolation<MasterPatch, SlavePatch>::"
            "rescaleWeightingFactors() const"
        )   << "Master or slave partially covered faces are not calculated."
            << abort(FatalError);
    }

    const labelList& partiallyUncoveredMasterFaces =
        *partiallyUncoveredMasterAddrPtr_;
    const labelList& partiallyUncoveredSlaveFaces =
        *partiallyUncoveredSlaveAddrPtr_;

    // Create a mask for partially covered master/slave faces
    boolList masterPCMask(maW.size(), false);
    boolList slavePCMask(saW.size(), false);

    forAll (partiallyUncoveredMasterFaces, pfmI)
    {
        masterPCMask[partiallyUncoveredMasterFaces[pfmI]] = true;
    }
    forAll (partiallyUncoveredSlaveFaces, pfsI)
    {
        slavePCMask[partiallyUncoveredSlaveFaces[pfsI]] = true;
    }

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
                "void GGIInterpolation<MasterPatch, SlavePatch>::"
                "rescaleWeightingFactors() const"
            )   << "Uncovered faces found.  On master: "
                << uncoveredMasterFaces().size()
                << " on slave: " << uncoveredSlaveFaces().size() << endl;
        }
      }

    forAll (saW, saWi)
    {
        scalar slaveWeightSum = Foam::sum(saW[saWi]);

        if (saW[saWi].size() > 0 && !slavePCMask[saWi])
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

        if (maW[maWi].size() > 0 && !masterPCMask[maWi])
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
GGIInterpolation<MasterPatch, SlavePatch>::findNonOverlappingFaces
(
    const scalarListList& patchWeights,
    const scalar& nonOverlapFaceTol   //  = min sum of the neighbour weights
) const
{
    tmp<labelField> tpatchFaceNonOverlapAddr(new labelField());
    labelField& patchFaceNonOverlapAddr = tpatchFaceNonOverlapAddr.ref();

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
        InfoIn("GGIInterpolation::findNonOverlappingFaces")
            << "   : Found " << patchFaceNonOverlapAddr.size()
            << " non-overlapping faces for this GGI patch" << endl;
    }

    return tpatchFaceNonOverlapAddr;
}


template<class MasterPatch, class SlavePatch>
void GGIInterpolation<MasterPatch, SlavePatch>::calcPartiallyCoveredFaces
(
    const scalarListList& patchWeights,
    const scalar& nonOverlapFaceTol,
    const bool isMaster
) const
{
    // Sanity checks first
    if
    (
        isMaster
     && (partiallyUncoveredMasterAddrPtr_ || masterFaceUncoveredFractionsPtr_)
    )
    {
        FatalErrorIn
        (
            "void GGIInterpolation<MasterPatch, SlavePatch>::"
            "calcPartiallyCoveredFaces() const"
        )   << "Already calculated master partially covered faces"
            << abort(FatalError);
    }

    if
    (
        !isMaster
     && (partiallyUncoveredSlaveAddrPtr_ || slaveFaceUncoveredFractionsPtr_)
    )
    {
        FatalErrorIn
        (
            "void GGIInterpolation<MasterPatch, SlavePatch>::"
            "calcPartiallyCoveredFaces() const"
        )   << "Already calculated slave partially covered faces"
            << abort(FatalError);
    }

    // Temporary storage
    DynamicList<label, 64> patchFacePartialOverlap(patchWeights.size());
    DynamicList<scalar, 64> uncoveredFaceFractions(patchWeights.size());

    // Scan the list of patch weights and collect ones inbetween
    // nonOverlapFaceTol and uncoveredFaceAreaTol_
    forAll (patchWeights, paWi)
    {
        const scalar sumWeightsFace = sum(patchWeights[paWi]);

        if
        (
            sumWeightsFace > nonOverlapFaceTol
         && sumWeightsFace <= uncoveredFaceAreaTol_
        )
        {
            // This is considered partially overlapped face, store the index and
            // the non-overlapping area (1 - sum of weights)
            patchFacePartialOverlap.append(paWi);
            uncoveredFaceFractions.append(1.0 - sumWeightsFace);
        }
    }

    // Transfer the storage
    if (isMaster)
    {
        // Allocate master side
        partiallyUncoveredMasterAddrPtr_ =
            new labelList(move(patchFacePartialOverlap));
        masterFaceUncoveredFractionsPtr_ =
            new scalarField(move(uncoveredFaceFractions));

        if (debug)
        {
            InfoIn("GGIInterpolation::calcPartiallyCoveredFaces")
                << "   : Found " << partiallyUncoveredMasterAddrPtr_->size()
                << " partially overlapping faces for master GGI patch" << endl;

            if (partiallyUncoveredMasterAddrPtr_->size())
            {
                Info<< "Max uncoverage: "
                    << max(*masterFaceUncoveredFractionsPtr_)
                    << ", min uncoverage: "
                    << min(*masterFaceUncoveredFractionsPtr_)
                    << endl;
            }
        }
    }
    else
    {
        // Allocate slave side
        partiallyUncoveredSlaveAddrPtr_ =
            new labelList(move(patchFacePartialOverlap));
        slaveFaceUncoveredFractionsPtr_ =
            new scalarField(move(uncoveredFaceFractions));

        if (debug)
        {
            InfoIn("GGIInterpolation::calcPartiallyCoveredFaces")
                << "   : Found " << partiallyUncoveredSlaveAddrPtr_->size()
                << " partially overlapping faces for slave GGI patch" << endl;

            if (partiallyUncoveredSlaveAddrPtr_->size())
            {
                Info<< "Max uncoverage: "
                    << max(*slaveFaceUncoveredFractionsPtr_)
                    << ", min uncoverage: "
                    << min(*slaveFaceUncoveredFractionsPtr_)
                    << endl;
            }
        }
    }
}


template<class FromPatch, class ToPatch>
void GGIInterpolation<FromPatch, ToPatch>::
calcMasterPointAddressing() const
{
    Info << "calcMasterPointAddressing() const" << endl;

    // Find master points addressing
    if (masterPointAddressingPtr_)
    {
        FatalErrorIn
        (
            "void ExtendedGGIInterpolation::"
            "calcMasterPointAddressing() const"
        )
            << "Master points addressing already exists"
                << abort(FatalError);
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

    const labelListList& masterFaceAddr = this->masterAddr();

    const labelListList& pointFaces = this->masterPatch().pointFaces();

    const faceList& slaveFaces = this->slavePatch().localFaces();
    const pointField& slavePoints = this->slavePatch().localPoints();

    const pointField& masterPoints = this->masterPatch().localPoints();

    forAll(masterPointAddr, pointI)
    {
        const point& P = masterPoints[pointI];

        labelHashSet possibleSlaveFacesSet;

        const labelList& curPointFaces = pointFaces[pointI];
        forAll(curPointFaces, faceI)
        {
            label curFace = curPointFaces[faceI];

            const labelList& curSlaveFaces = masterFaceAddr[curFace];

            forAll(curSlaveFaces, fI)
            {
                if (!possibleSlaveFacesSet.found(curSlaveFaces[fI]))
                {
                    possibleSlaveFacesSet.insert(curSlaveFaces[fI]);
                }
            }
        }

        labelList possibleSlaveFaces = possibleSlaveFacesSet.toc();

        scalar MinEta = -GREAT;
        labelPair faceTriangle(-1, -1);
        scalar distance = GREAT;

        forAll(possibleSlaveFaces, faceI)
        {
            label curSlaveFace = possibleSlaveFaces[faceI];

            const face& f = slaveFaces[curSlaveFace];

            point ctr = Foam::average(f.points(slavePoints));

            point nextPoint = ctr;

            for (label pI = 0; pI < f.size(); pI++)
            {
                nextPoint = slavePoints[f.nextLabel(pI)];

                triPointShapeRef t
                (
                    slavePoints[f[pI]],
                    nextPoint,
                    ctr
                );

                vector n = t.normal();
                scalar A = mag(n);
                n /= A;

                // Intersection point
                point I = P + n*(n&(t.a() - P));

                // Areal coordinates
                scalarField eta(3, 0);

                eta[0] = (triPointRef(I, t.b(), t.c()).normal() & n)/A;
                eta[1] = (triPointRef(I, t.c(), t.a()).normal() & n)/A;
                eta[2] = (triPointRef(I, t.a(), t.b()).normal() & n)/A;

                scalar minEta = min(eta);

                if (minEta > MinEta)
                {
                    MinEta = minEta;
                    faceTriangle.first() = curSlaveFace;
                    faceTriangle.second() = pI;

                    distance = ((P - I)&n);
                }
            }
        }

        masterPointAddr[pointI] = faceTriangle;
        masterPointDist[pointI] = distance;
    }

    Info << "Extended GGI, master point distance, max: "
        << max(masterPointDist)
        << ", avg: " << average(masterPointDist)
        << ", min: " << min(masterPointDist) << endl;


    // Check orientation

    const pointField& masterPointNormals =
        this->masterPatch().pointNormals();

    const vectorField& slaveFaceNormals =
        this->slavePatch().faceNormals();

    scalarField orientation(masterPointAddr.size(), 0);

    label nIncorrectPoints = 0;

    forAll(masterPointAddr, pointI)
    {
        orientation[pointI] =
            (
                masterPointNormals[pointI]
              & slaveFaceNormals[masterPointAddr[pointI].first()]
            );

        if (orientation[pointI] > -SMALL)
        {
            nIncorrectPoints++;
        }
    }

//     Info << "Extended GGI, master point orientation (<0), max: "
//         << max(orientation)
//         << ", min: " << min(orientation) << ", nIncorrectPoints: "
//         << nIncorrectPoints << "/" << masterPointAddr.size() << endl;
}

template<class FromPatch, class ToPatch>
void GGIInterpolation<FromPatch, ToPatch>::
calcMasterPointWeights() const
{
    // Find master point weights
    if (masterPointWeightsPtr_)
    {
        FatalErrorIn
        (
            "void ExtendedGGIInterpolation::"
            "calcMasterPointAddressing() const"
        )
            << "Master point weights already exist"
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

            triPointShapeRef t
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
void GGIInterpolation<FromPatch, ToPatch>::
calcSlavePointAddressing() const
{
    Info << "calcSlavePointAddressing() const" << endl;

    // Find master points addressing
    if (slavePointAddressingPtr_)
    {
        FatalErrorIn
        (
            "void ExtendedGGIInterpolation::"
            "calcSlavePointAddressing() const"
        )
            << "Slave points addressing already exists"
                << abort(FatalError);
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

    const labelListList& slaveFaceAddr = this->slaveAddr();

    const labelListList& pointFaces = this->slavePatch().pointFaces();

    const faceList& masterFaces = this->masterPatch().localFaces();
    const pointField& masterPoints = this->masterPatch().localPoints();

    const pointField& slavePoints = this->slavePatch().localPoints();

    forAll(slavePointAddr, pointI)
    {
        const point& P = slavePoints[pointI];

        labelHashSet possibleMasterFacesSet;

        const labelList& curPointFaces = pointFaces[pointI];
        forAll(curPointFaces, faceI)
        {
            label curFace = curPointFaces[faceI];

            const labelList& curMasterFaces = slaveFaceAddr[curFace];

            forAll(curMasterFaces, fI)
            {
                if (!possibleMasterFacesSet.found(curMasterFaces[fI]))
                {
                    possibleMasterFacesSet.insert(curMasterFaces[fI]);
                }
            }
        }

        labelList possibleMasterFaces = possibleMasterFacesSet.toc();

        scalar MinEta = -GREAT;
        labelPair faceTriangle(-1, -1);
        scalar distance = GREAT;

        forAll(possibleMasterFaces, faceI)
        {
            label curMasterFace = possibleMasterFaces[faceI];

            const face& f = masterFaces[curMasterFace];

            point ctr = Foam::average(f.points(masterPoints));

            point nextPoint = ctr;

            for (label pI = 0; pI < f.size(); pI++)
            {
                nextPoint = masterPoints[f.nextLabel(pI)];

                triPointRef t
                (
                    masterPoints[f[pI]],
                    nextPoint,
                    ctr
                );

                vector n = t.normal();
                scalar A = mag(n);
                n /= A;

                // Intersection point
                point I = P + n*(n&(t.a() - P));

                // Areal coordinates
                scalarField eta(3, 0);

                eta[0] = (triPointRef(I, t.b(), t.c()).normal() & n)/A;
                eta[1] = (triPointRef(I, t.c(), t.a()).normal() & n)/A;
                eta[2] = (triPointRef(I, t.a(), t.b()).normal() & n)/A;

                scalar minEta = min(eta);

                if (minEta > MinEta)
                {
                    MinEta = minEta;
                    faceTriangle.first() = curMasterFace;
                    faceTriangle.second() = pI;

                    distance = ((P - I)&n);
                }
            }
        }

//         Info << "MinEta " << MinEta << endl;

        slavePointAddr[pointI] = faceTriangle;
        slavePointDist[pointI] = distance;

//         Info << "slave " << pointI << ", "
//             << slavePointAddr[pointI] << endl;
    }

    Info << "Extended GGI, slave point distance, max: "
        << max(slavePointDist)
        << ", avg: " << average(slavePointDist)
        << ", min: " << min(slavePointDist) << endl;

    // Check orientation

    const pointField& slavePointNormals =
        this->slavePatch().pointNormals();

    const vectorField& masterFaceNormals =
        this->masterPatch().faceNormals();

    scalarField orientation(slavePointAddr.size(), 0);

    label nIncorrectPoints = 0;

    forAll(slavePointAddr, pointI)
    {
        if (slavePointAddr[pointI].first() != -1)
        {
            orientation[pointI] =
                (
                    slavePointNormals[pointI]
                  & masterFaceNormals[slavePointAddr[pointI].first()]
                );

            if (orientation[pointI] > -SMALL)
            {
                nIncorrectPoints++;
            }
        }
    }

    Info << "Extended GGI, slave point orientation (<0), max: "
        << max(orientation)
        << ", min: " << min(orientation) << ", nIncorrectPoints: "
        << nIncorrectPoints << "/" << slavePointAddr.size() << endl;
}

template<class FromPatch, class ToPatch>
void GGIInterpolation<FromPatch, ToPatch>::
calcSlavePointWeights() const
{
    // Find master point weights
    if (slavePointWeightsPtr_)
    {
        FatalErrorIn
        (
            "void ExtendedGGIInterpolation::"
            "calcSlavePointAddressing() const"
        )
            << "Slave point weights already exist"
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

            triPointShapeRef t
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



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
