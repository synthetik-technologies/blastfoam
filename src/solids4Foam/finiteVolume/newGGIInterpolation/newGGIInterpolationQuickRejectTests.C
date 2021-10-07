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
    3D and 2D quick reject tests for filtering out the neighborhood
    of the GGI master patch faces

Author
    Martin Beaudoin, Hydro-Quebec, (2008)
    updateNeighboursAABB added by
    Philip Cardiff, UCD
    Tian Tang, Bekaert
    Peter De Jaeger, Bekaert

\*---------------------------------------------------------------------------*/

#include "boundBox.H"
#include "plane.H"
#include "transformField.H"
#include "octree.H"
#include "octreeDataBoundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
const Foam::debug::tolerancesSwitch
newGGIInterpolation<MasterPatch, SlavePatch>::faceBoundBoxExtendSpanFraction_
(
    "GGIFaceBoundBoxExtendSpanFraction",
    1.0e-2,
    "GGI faces bounding box expansion factor. "
    "Add robustness for quick-search algo. Keep it to a few percent."
);


template<class MasterPatch, class SlavePatch>
const Foam::debug::optimisationSwitch
newGGIInterpolation<MasterPatch, SlavePatch>::octreeSearchMinNLevel_
(
    "GGIOctreeSearchMinNLevel",
    3,
    "GGI neighbouring facets octree-based search: "
    "minNlevel parameter for octree"
);


template<class MasterPatch, class SlavePatch>
const Foam::debug::optimisationSwitch
newGGIInterpolation<MasterPatch, SlavePatch>::octreeSearchMaxLeafRatio_
(
    "GGIOctreeSearchMaxLeafRatio",
    3,
    "GGI neighbouring facets octree-based search: "
    "maxLeafRatio parameter for octree"
);


template<class MasterPatch, class SlavePatch>
const Foam::debug::optimisationSwitch
newGGIInterpolation<MasterPatch, SlavePatch>::octreeSearchMaxShapeRatio_
(
    "GGIOctreeSearchMaxShapeRatio",
    1,
    "GGI neighbouring facets octree-based search: "
    "maxShapeRatio parameter for octree"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// From: http://www.gamasutra.com/features/20000330/bobic_02.htm
// One of the most primitive ways of doing collision detection is to
// approximate each object or a part of the object with a sphere, and
// then check whether spheres intersect each other. This method is
// widely used even today because it's computationally inexpensive.
// We merely check whether the distance between the centers of two
// spheres is less than the sum of the two radii (which indicates that
// a collision has occurred). Even better, if we calculate whether the
// distance squared is less than the sum of the radii squared, then we
// eliminate that nasty square root in our distance calculation, but
// only if we take care of this:
//
// In the interval [0,1], sqr(x) underestimate x, so we will have
// false negative for the neighbors, this is bad...  In the interval
// ]1, infinite], sqr(x) will overestimate x, so we will get false
// positive, which is not that bad for a quick reject test
//
// So we will check the values of the sqr(x) before doing comparisons,
// if sqr(x) is < 1, then we will pay for the sqrt. If sqr(x) > 1,
// then we only compare squared values.
//
// Of course, we will end up comparing squared and and no-squared
// values, which is a bit weird; but we must remember that this is a
// quick reject test, and nothing more.
//
// Overall, we will overestimate the number of neighbours, which is
// much more preferable than to underestimate. We will have to use
// another test to remove the false positives (Separating Axis Theorem
// test)
//
// This constitutes the first good quick reject test before going into
// more computing intensive tasks with the calculation of the GGI
// weights

template<class MasterPatch, class SlavePatch>
void newGGIInterpolation<MasterPatch, SlavePatch>::findNeighbours3D
(
    labelListList& result
) const
{
    // Allocation to local size.  HJ, 27/Apr/2016
    List<DynamicList<label, 8> > candidateMasterNeighbors(parMasterSize());

    // First, compute the face center and the sphere radius (squared)
    // of the slave patch faces so we will not have to recompute this
    // n times
    vectorField slaveFaceCentre(slavePatch_.size());
    scalarField slaveRadius2(slavePatch_.size());

    // Since the parallelised search is done on the master, all of
    // slave face data is needed.  HJ, 27/Apr/2016
    forAll (slavePatch_, faceSi)
    {
        // Grab points from slave face
        pointField curFacePoints =
            slavePatch_[faceSi].points(slavePatch_.points());

        slaveFaceCentre[faceSi] =
            slavePatch_[faceSi].centre(slavePatch_.points());

        if (doTransform())
        {
            // Transform points and normal to master plane

            if (forwardT_.size() == 1)
            {
                // Constant transform
                transform
                (
                    curFacePoints,
                    forwardT_[0],
                    curFacePoints
                );

                slaveFaceCentre[faceSi] = transform
                (
                    forwardT_[0],
                    slaveFaceCentre[faceSi]
                );
            }
            else
            {
                transform
                (
                    curFacePoints,
                    forwardT_[faceSi],
                    curFacePoints
                );

                slaveFaceCentre[faceSi] = transform
                (
                    forwardT_[faceSi],
                    slaveFaceCentre[faceSi]
                );
            }
        }

        boundBox bbSlave(curFacePoints, false);

        scalar tmpValue = 0.25*Foam::magSqr(bbSlave.max() - bbSlave.min());

        // We compare squared distances, so save the sqrt() if value > 1.
        if (tmpValue < 1)
        {
            // Take the sqrt, otherwise, we underestimate the radius
            slaveRadius2[faceSi] = sqrt(tmpValue);
        }
        else
        {
            slaveRadius2[faceSi] = tmpValue;
        }
    }

    // Next, we search for each master face a list of potential neighbors

    // Parallel search split.  HJ, 27/Apr/2016
    const label pmStart = this->parMasterStart();

    for
    (
        label faceMi = pmStart;
        faceMi < this->parMasterEnd();
        faceMi++
    )
//     forAll (masterPatch_, faceMi)
    {
        // For each masterPatch faces, compute the bounding box. With
        // this, we compute the radius of the bounding sphere for this
        // face
        boundBox bbMaster
        (
            masterPatch_[faceMi].points(masterPatch_.points()),
            false
        );

        // We will compare squared distances, so save the sqrt() if > 1.0.
        scalar masterRadius2 =
            Foam::magSqr(bbMaster.max() - bbMaster.min())/4.0;

        // Take the sqrt, otherwise, we underestimate the radius
        if (masterRadius2 < 1.0)
        {
            masterRadius2 = sqrt(masterRadius2);
        }

        // This is the inner loop, so the face centres and face radii
        // for the slave faces are already pre-computed
        forAll (slavePatch_, faceSi)
        {
            // We test if the squared distance between the two face
            // centers is less than the sum of the 2 squared radii from
            // each face bounding sphere.

            // If so, we have found possibly 2 neighbours
            scalar distFaceCenters =
                Foam::magSqr
                (
                    masterPatch_[faceMi].centre(masterPatch_.points())
                  - slaveFaceCentre[faceSi]
                );

            if (distFaceCenters < 1.0)
            {
                distFaceCenters = sqrt(distFaceCenters);
            }

            if (distFaceCenters < (masterRadius2 + slaveRadius2[faceSi]))
            {
                candidateMasterNeighbors[faceMi - pmStart].append(faceSi);
            }
        }
    }

    // Repack the list.  Local size
    result.setSize(parMasterSize());

    // Parallel search split: local size.  HJ, 27/Apr/2016
    forAll (result, i)
    {
        result[i].transfer(candidateMasterNeighbors[i].shrink());
    }

    // Note: since the neighbours are used to perform a search, there is
    // no need to do a global reduce of the candidates.  The rest of the list
    // (searched on other processors) will remain unused.
    // HJ, 27/Apr/2016
}


// This algorithm find the faces in proximity of another face based
// on the AABB (Axis Aligned Boundary Box. These tests are done in 3D
// space.
//
// For a given face and a potential neighbours, we construct first an
// "augmented BB" by substracting/adding the slave delta BB to the
// BBmin/BBmax of the master face. Then, we compare if the "augmented
// BB" intersects with the potential neighbour BB.
//
// This is a fairly quick reject test, based only on additions and
// comparisons... But still a n^2 test... very costly
//
// Furthermore, before selecting a given face as a potential
// neighbour, we evaluate the featureCos of both master and slave
// face normals. We reject immediatly all the faces that might fall
// into the bounding box, but where the normals are too dissimilar.

template<class MasterPatch, class SlavePatch>
void newGGIInterpolation<MasterPatch, SlavePatch>::findNeighboursAABB
(
    labelListList& result
) const
{
    // Allocation to local size.  HJ, 27/Apr/2016
    List<DynamicList<label, 8> > candidateMasterNeighbors(parMasterSize());

    // Allocation to local size.  HJ, 27/Apr/2016
    List<boundBox> masterPatchBB(parMasterSize());

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
        masterPatchBB[faceMi - pmStart] = boundBox
        (
            masterPatch_[faceMi].points(masterPatch_.points()),
            false
        );
    }

    // Grab the slave patch faces bounding boxes, plus compute its
    // extend or delta
    List<boundBox> slavePatchBB(slavePatch_.size());
    pointField deltaBBSlave(slavePatch_.size());

    // We expect that any neighbour face to face intersection will fall
    // within augmented BB.
    vectorField slaveFaceBBminThickness(slavePatch_.size(), vector::zero);

    const faceList& slaveLocalFaces = slavePatch_.localFaces();
    vectorField slaveNormals = slavePatch_.faceNormals();
    const pointField& slaveLocalPoints = slavePatch_.localPoints();

    // Transform slave normals to master plane if needed
    if (doTransform())
    {
        if (forwardT_.size() == 1)
        {
            transform(slaveNormals, forwardT_[0], slaveNormals);
        }
        else
        {
            transform(slaveNormals, forwardT_, slaveNormals);
        }
    }

    forAll (slaveFaceBBminThickness, sI)
    {
        scalar maxEdgeLength = 0.0;

        // Let's use the length of the longest edge from each faces
        edgeList el = slaveLocalFaces[sI].edges();

        forAll (el, elI)
        {
            scalar edgeLength = el[elI].mag(slaveLocalPoints);
            maxEdgeLength = Foam::max(edgeLength, maxEdgeLength);
        }

        // Make sure our offset is positive. Ugly, but cheap
        vector posNormal = cmptMag(slaveNormals[sI]);

        slaveFaceBBminThickness[sI] = posNormal*maxEdgeLength;
    }

    // Iterate over slave patch faces, compute its bounding box,
    // using a possible transformation and separation for cyclic patches
    forAll (slavePatch_, faceSi)
    {
        pointField curFacePoints =
            slavePatch_[faceSi].points(slavePatch_.points());

        if (doTransform())
        {
            if (forwardT_.size() == 1)
            {
                transform(curFacePoints, forwardT_[0], curFacePoints);
            }
            else
            {
                transform(curFacePoints, forwardT_[faceSi], curFacePoints);
            }
        }

        if (doSeparation())
        {
            if (forwardSep_.size() == 1)
            {
                curFacePoints += forwardSep_[0];
            }
            else
            {
                curFacePoints += forwardSep_[faceSi];
            }
        }

        slavePatchBB[faceSi] = boundBox(curFacePoints, false);

        // We compute the extent of the slave face BB.
        // Plus, we boost it a little bit, just to stay clear
        // of floating point numerical issues when doing intersections
        // Let's boost by 10%.
        deltaBBSlave[faceSi] =
            1.1*
            (
                slavePatchBB[faceSi].max()
              - slavePatchBB[faceSi].min()
              + slaveFaceBBminThickness[faceSi]
            );
    }

    // Visit each master patch face BB,
    // augment it with the info from the slave patch face BB
    // then, compute the intersection
    const vectorField& masterFaceNormals = masterPatch_.faceNormals();

    // Check which faces are in the region of interest
    boolList checkMasterFace(masterPatchBB.size(), false);
    forAll (masterPatchBB, faceMi)
    {
        if (regionOfInterest_.contains(masterPatchBB[faceMi].midpoint()))
        {
            checkMasterFace[faceMi] = true;
        }
    }

    boolList checkSlaveFace(slavePatchBB.size(), false);
    forAll (slavePatchBB, faceSi)
    {
        if (regionOfInterest_.contains(slavePatchBB[faceSi].midpoint()))
        {
            checkSlaveFace[faceSi] = true;
        }
    }

    int countMaster(0);
    int countSlave(0);

    // Parallel search split.  HJ, 27/Apr/2016
    for
    (
        label faceMi = this->parMasterStart();
        faceMi < this->parMasterEnd();
        faceMi++
    )
    //forAll (masterPatch_, faceMi)
    {
        if (checkMasterFace[faceMi - pmStart])
        {
            countMaster++;
            countSlave = 0;

            forAll(slavePatchBB, faceSi)
            {
                if (checkSlaveFace[faceSi])
                {
                    countSlave++;

                    // Compute the augmented AABB
                    boundBox augmentedBBMaster
                    (
                        masterPatchBB[faceMi - pmStart].min()
                      - deltaBBSlave[faceSi],
                        masterPatchBB[faceMi - pmStart].max()
                      + deltaBBSlave[faceSi]
                    );

                    if (augmentedBBMaster.overlaps(slavePatchBB[faceSi]))
                    {
                        // Compute featureCos between the two face normals
                        // before adding to the list of candidates
                        scalar featureCos =
                            masterFaceNormals[faceMi] & slaveNormals[faceSi];

                        if (mag(featureCos) > featureCosTol_)
                        {
                            candidateMasterNeighbors[faceMi - pmStart].append
                            (
                                faceSi
                            );
                        }
                    }
                }
            }
        }
    }

    // Repack the list.  Local size
    result.setSize(parMasterSize());

    // Parallel search split: local size.  HJ, 27/Apr/2016
    forAll (result, i)
    {
        result[i].transfer(candidateMasterNeighbors[i].shrink());
    }
}


// This algorithm find the faces in proximity of another face based
// on the face BB (Bounding Box) and an octree of bounding boxes.
template<class MasterPatch, class SlavePatch>
void newGGIInterpolation<MasterPatch, SlavePatch>::findNeighboursBBOctree
(
    labelListList& result
) const
{
    // Allocation to local size.  HJ, 27/Apr/2016
    List<DynamicList<label, 8> > candidateMasterNeighbors(parMasterSize());

    // Allocation to local size.  HJ, 27/Apr/2016
    treeBoundBoxList lmasterFaceBB(parMasterSize());

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
        pointField facePoints
        (
            masterPatch_[faceMi].points(masterPatch_.points())
        );

        // Construct face BB with an extension of face span defined by the
        //  global tolerance factor faceBoundBoxExtendSpanFraction_
        // (1% by default)
        treeBoundBox bbFaceMaster(facePoints);

        lmasterFaceBB[faceMi - pmStart] =
            bbFaceMaster.extend(faceBoundBoxExtendSpanFraction_());
    }

    // Initialize the list of slave patch faces bounding box
    treeBoundBoxList lslaveFaceBB(slavePatch_.size());

    forAll (slavePatch_, faceSi)
    {
        pointField facePoints
        (
            slavePatch_[faceSi].points(slavePatch_.points())
        );

        // possible transformation and separation for cyclic patches
        if (doTransform())
        {
            if (forwardT_.size() == 1)
            {
                transform(facePoints, forwardT_[0], facePoints);
            }
            else
            {
                transform(facePoints, forwardT_[faceSi], facePoints);
            }
        }

        if (doSeparation())
        {
            if (forwardSep_.size() == 1)
            {
                facePoints += forwardSep_[0];
            }
            else
            {
                facePoints += forwardSep_[faceSi];
            }
        }

        // Construct face BB with an extension of face span defined by the
        //  global tolerance factor faceBoundBoxExtendSpanFraction_
        // (1% by default)
        treeBoundBox bbFaceSlave(facePoints);

        lslaveFaceBB[faceSi] =
            bbFaceSlave.extend(faceBoundBoxExtendSpanFraction_());
    }

    // Create the slave octreeData, using the boundBox flavor
    octreeDataBoundBox slaveDataBB(lslaveFaceBB);

    // Overall slave patch BB
    treeBoundBox slaveOverallBB(slavePatch_.points());

    // Create the slave patch octree

    octree<octreeDataBoundBox> slavePatchOctree
    (
        slaveOverallBB,              // overall search domain
        slaveDataBB,
        octreeSearchMinNLevel_(),      // min number of levels
        octreeSearchMaxLeafRatio_(),   // max avg. size of leaves
        octreeSearchMaxShapeRatio_()   // max avg. duplicity.
    );

    const vectorField& masterFaceNormals = masterPatch_.faceNormals();
    vectorField slaveNormals = slavePatch_.faceNormals();

    // Transform slave normals to master plane if needed
    if (doTransform())
    {
        if (forwardT_.size() == 1)
        {
            transform(slaveNormals, forwardT_[0], slaveNormals);
        }
        else
        {
            transform(slaveNormals, forwardT_, slaveNormals);
        }
    }

    // Visit each master patch face BB and find the potential neighbours
    // using the slave patch octree

    // Parallel search split.  HJ, 27/Apr/2016
    for
    (
        label faceMi = this->parMasterStart();
        faceMi < this->parMasterEnd();
        faceMi++
    )
//     forAll (masterPatch_, faceMi)
    {
        // List of candidate neighbours
        labelList overlappedFaces  =
            slavePatchOctree.findBox(lmasterFaceBB[faceMi - pmStart]);

        forAll (overlappedFaces, ovFi)
        {
            label faceSi = overlappedFaces[ovFi];

            // Compute and verify featureCos between the two face normals
            // before adding to the list of candidates
            scalar featureCos =
                masterFaceNormals[faceMi] & slaveNormals[faceSi];

            if (mag(featureCos) > featureCosTol_)
            {
                candidateMasterNeighbors[faceMi - pmStart].append(faceSi);
            }
        }
    }

    // Repack the list.  Local size
    result.setSize(parMasterSize());

    // Parallel search split: local size.  HJ, 27/Apr/2016
    forAll (result, i)
    {
        result[i].transfer(candidateMasterNeighbors[i].shrink());
    }

    // Note: since the neighbours are used to perform a search, there is
    // no need to do a global reduce of the candidates.  The rest of the list
    // (searched on other processors) will remain unused.
    // HJ, 27/Apr/2016
}


template<class MasterPatch, class SlavePatch>
void newGGIInterpolation<MasterPatch, SlavePatch>::updateNeighboursAABB
(
    labelListList& result
) const
{
    // Allocation to local size.  HJ, 27/Apr/2016
    List<DynamicList<label, 8> > candidateMasterNeighbors(parMasterSize());


    // Create master face bounding boxes


    // Allocation to local size.  HJ, 27/Apr/2016
    List<boundBox> masterPatchBB(parMasterSize());

    // Parallel search split.  HJ, 27/Apr/2016
    const label pmStart = parMasterStart();
    const label pmEnd = parMasterEnd();

    for (label faceMi = pmStart; faceMi < pmEnd; faceMi++)
    {
        masterPatchBB[faceMi - pmStart] = boundBox
        (
            masterPatch_[faceMi].points(masterPatch_.points()),
            false
        );
    }


    // Create augmented slave face bounding boxes


    // Grab the slave patch faces bounding boxes, plus compute its
    // extend or delta
    List<boundBox> slavePatchBB(slavePatch_.size());
    pointField deltaBBSlave(slavePatch_.size());

    // We expect that any neighbour face to face intersection will fall
    // within augmented BB.
    vectorField slaveFaceBBminThickness(slavePatch_.size(), vector::zero);

    const faceList& slaveLocalFaces = slavePatch_.localFaces();
    vectorField slaveNormals = slavePatch_.faceNormals();
    const pointField& slaveLocalPoints = slavePatch_.localPoints();

    // Transform slave normals to master plane if needed
    if (doTransform())
    {
        if (forwardT_.size() == 1)
        {
            transform(slaveNormals, forwardT_[0], slaveNormals);
        }
        else
        {
            transform(slaveNormals, forwardT_, slaveNormals);
        }
    }

    forAll(slaveFaceBBminThickness, sI)
    {
        scalar maxEdgeLength = 0.0;

        // Let's use the length of the longest edge from each faces
        const edgeList el = slaveLocalFaces[sI].edges();

        forAll(el, elI)
        {
            const scalar edgeLength = el[elI].mag(slaveLocalPoints);
            maxEdgeLength = Foam::max(edgeLength, maxEdgeLength);
        }

        // Make sure our offset is positive. Ugly, but cheap
        const vector posNormal = cmptMag(slaveNormals[sI]);

        slaveFaceBBminThickness[sI] = posNormal*maxEdgeLength;
    }

    // Iterate over slave patch faces, compute its bounding box,
    // using a possible transformation and separation for cyclic patches
    forAll(slavePatch_, faceSi)
    {
        pointField curFacePoints =
            slavePatch_[faceSi].points(slavePatch_.points());

        if (doTransform())
        {
            if (forwardT_.size() == 1)
            {
                transform(curFacePoints, forwardT_[0], curFacePoints);
            }
            else
            {
                transform(curFacePoints, forwardT_[faceSi], curFacePoints);
            }
        }

        if (doSeparation())
        {
            if (forwardSep_.size() == 1)
            {
                curFacePoints += forwardSep_[0];
            }
            else
            {
                curFacePoints += forwardSep_[faceSi];
            }
        }

        slavePatchBB[faceSi] = boundBox(curFacePoints, false);

        // We compute the extent of the slave face BB.
        // Plus, we boost it a little bit, just to stay clear
        // of floating point numerical issues when doing intersections
        // Let's boost by 10%.
        // PC: we'll boost by a larger factor as faces could be in contact by a
        // large ammount
        deltaBBSlave[faceSi] =
            1.1*
            //2.0*
            (
                slavePatchBB[faceSi].max()
              - slavePatchBB[faceSi].min()
              + slaveFaceBBminThickness[faceSi]
            );
    }


    // Method: Philip Cardiff and Tian Tang
    // Starting from the previous neighbours, we will perform a walk from
    // face-to-face to find new potential neighbours: this method takes N time.


    const vectorField& masterFaceNormals = masterPatch_.faceNormals();
    const labelListList& masterFaceFaces = masterPatch_.faceFaces();
    const labelListList& slaveFaceFaces = slavePatch_.faceFaces();

    // Check if there are any neighbours stored from the previous iteration
    if (prevCandidateMasterNeighbors_.size() == 0)
    {
        Info<< "    " << typeName << " : n2 search" << endl;

        // Set the size of the stored neighbours list
        prevCandidateMasterNeighbors_.setSize(parMasterSize());

        for (label faceMi = pmStart; faceMi < pmEnd; faceMi++)
        {
            forAll(slavePatchBB, faceSi)
            {
                // Compute the augmented AABB
                boundBox augmentedBBMaster
                (
                    masterPatchBB[faceMi - pmStart].min()
                  - deltaBBSlave[faceSi],
                    masterPatchBB[faceMi - pmStart].max()
                  + deltaBBSlave[faceSi]
                );

                if (augmentedBBMaster.overlaps(slavePatchBB[faceSi]))
                {
                    // Compute featureCos between the two face normals
                    // before adding to the list of candidates
                    const scalar featureCos =
                        masterFaceNormals[faceMi] & slaveNormals[faceSi];

                    if (mag(featureCos) > featureCosTol_)
                    {
                        candidateMasterNeighbors[faceMi - pmStart].append
                        (
                            faceSi
                        );
                    }
                }
            }
        }
    }
    else
    {
        // We will use the previous neighbours as a starting guess to find the
        // new candidate neighbours
        for (label faceMi = pmStart; faceMi < pmEnd; faceMi++)
        {
            // Neighbour for this master face from the last iteration
            // Note: this array becomes invalid below
            labelList& curPrevResult =
                prevCandidateMasterNeighbors_[faceMi - pmStart];

            // If this master face has no previous neighbours then we
            // will use the previous neighbours of the face-faces as a guess
            if (curPrevResult.size() == 0)
            {
                labelHashSet prevNeiSet;
                const labelList& curFaceFaces = masterFaceFaces[faceMi];
                forAll(curFaceFaces, ffI)
                {
                    const label curFaceID = curFaceFaces[ffI];

                    if (curFaceID >= pmStart && curFaceID < pmEnd)
                    {
                        if
                        (
                            prevCandidateMasterNeighbors_
                            [
                                curFaceID - pmStart
                            ].size() > 0
                        )
                        {
                            if (!prevNeiSet.found(curFaceID))
                            {
                                prevNeiSet.insert(curFaceID);
                            }
                        }
                    }
                }

                // Update the previous neighbours for this face
                curPrevResult = prevNeiSet.toc();
            }

            if (curPrevResult.size() > 0)
            {
                // Keep track of the slave faces that have been checked
                labelHashSet slaveFacesChecked;
                forAll(curPrevResult, fI)
                {
                    slaveFacesChecked.insert(curPrevResult[fI]);
                }

                // Transfer curPrevResult to a dynamic list
                DynamicList<label> facesToCheck(10*curPrevResult.size());
                facesToCheck.transfer(curPrevResult);

                do
                {
                    // Get the next face to check and remove it from the list
                    const label faceSi = facesToCheck.remove();

                    // Mark this face as having been checked
                    slaveFacesChecked.insert(faceSi);

                    // Compute the augmented AABB for the master face
                    boundBox augmentedBBMaster
                    (
                        masterPatchBB[faceMi - pmStart].min()
                      - deltaBBSlave[faceSi],
                        masterPatchBB[faceMi - pmStart].max()
                      + deltaBBSlave[faceSi]
                    );

                    // Check if the master and slave BB overlap
                    if (augmentedBBMaster.overlaps(slavePatchBB[faceSi]))
                    {
                        // Compute featureCos between the two face normals
                        // before adding to the list of candidates
                        const scalar featureCos =
                            masterFaceNormals[faceMi] & slaveNormals[faceSi];

                        if (mag(featureCos) > featureCosTol_)
                        {
                            // Add the face to the list of candidate neighbours
                            candidateMasterNeighbors[faceMi - pmStart].append
                            (
                                faceSi
                            );

                            // Add face-face neighbours to be checked
                            const labelList& curSlaveFaceFaces =
                                slaveFaceFaces[faceSi];

                            forAll(curSlaveFaceFaces, ffI)
                            {
                                const label faceID = curSlaveFaceFaces[ffI];

                                if (!slaveFacesChecked.found(faceID))
                                {
                                    facesToCheck.append(faceID);
                                }
                            }
                        }
                    }
                } while (facesToCheck.size());
            }
        }
    }

    // Repack the list.  Local size
    result.setSize(parMasterSize());

    // Parallel search split: local size.  HJ, 27/Apr/2016
    forAll(result, i)
    {
        result[i].transfer(candidateMasterNeighbors[i].shrink());
    }

    // Update the previous neighbours to be used for the search next time
    prevCandidateMasterNeighbors_ = result;
}


// Projects a list of points onto a plane located at planeOrig,
// oriented along planeNormal.  Return the projected points in a
// pointField, and the normal distance of each points from the
// projection plane
template<class MasterPatch, class SlavePatch>
tmp<pointField>
newGGIInterpolation<MasterPatch, SlavePatch>::projectPointsOnPlane
(
    const pointField& lpoints,
    const vector& planeOrig,
    const vector& planeDirection,
    scalarField& distanceProjection
) const
{
    tmp<pointField> tprojectedPoints(new pointField(lpoints.size()));
    pointField& projectedPoints = tprojectedPoints();

    vector normalVector = planeDirection/(mag(planeDirection) + VSMALL);

    scalarField dist(lpoints.size(), 0.0);

    if (lpoints.size() > 3 && mag(normalVector) > SMALL)
    {
        // Construct the plane
        plane projectionPlane(planeOrig, normalVector);

        forAll (lpoints, pointI)
        {
            projectedPoints[pointI] =
                projectionPlane.nearestPoint(lpoints[pointI]);

            dist[pointI] = projectionPlane.distance(lpoints[pointI]);
        }
    }
    else
    {
          // Triangle... nothing to project, the points are already
          // located on the right plane
        projectedPoints = lpoints;
    }

    // Transfer the projection distance
    distanceProjection = dist;

    return tprojectedPoints;
}


// Compute an orthonormal basis (u, v, w) where: w is aligned on the
// normalVector u is the direction from normalVectorCentre to the most
// distant point in the list pointsOnPlane v = w^u
//
// Everything is normalized
template<class MasterPatch, class SlavePatch>
typename newGGIInterpolation<MasterPatch, SlavePatch>::orthoNormalBasis
newGGIInterpolation<MasterPatch, SlavePatch>::computeOrthonormalBasis
(
    const vector& normalVectorCentre,
    const vector& normalVector,
    const pointField& pointsOnPlane
) const
{
    // The orthonormal basis uvw
    orthoNormalBasis uvw;

    vector u = pTraits<vector>::zero;
    vector v = pTraits<vector>::zero;
    vector w = normalVector;

    // Normalized w
    w /= mag(w) + VSMALL;

    // Find the projected point from the master face that is the most distant
    scalar longestDistanceFromCenter = 0;

    forAll (pointsOnPlane, pointI)
    {
        vector delta = pointsOnPlane[pointI] - normalVectorCentre;

        if (mag(delta) > longestDistanceFromCenter)
        {
            longestDistanceFromCenter = mag(delta);
            u = delta; // This is our next candidate for the u direction
        }
    }

    // Normalized u
    u /= mag(u) + VSMALL;

    // Compute v from u and w
    v = w ^ u;    // v = w^u;

    // We got ourselves a new orthonormal basis to play with
    uvw[0] = u;
    uvw[1] = v;
    uvw[2] = w;

    return uvw;
}


// Projection of 3D points onto a 2D UV plane defined by an
// orthonormal basis We project onto the uv plane. Could be
// parametrized if we ever need to project onto other uvw planes
template<class MasterPatch, class SlavePatch>
List<point2D> newGGIInterpolation<MasterPatch, SlavePatch>::projectPoints3Dto2D
(
    const orthoNormalBasis& orthoBase,
    const vector& orthoBaseOffset,
    const pointField& pointsIn3D,
    scalarField& distanceProjection
) const
{
    List<point2D> pointsIn2D(pointsIn3D.size());
    scalarField  dist(pointsIn3D.size(), 0.0);

    pointField pointsIn3DTranslated = pointsIn3D - orthoBaseOffset;

    // Project onto the uv plane. The distance from the plane is
    // computed with w
    forAll (pointsIn3D, pointsI)
    {
        // u component
        pointsIn2D[pointsI][0] =
            pointsIn3DTranslated[pointsI] & orthoBase[0];

        // v component
        pointsIn2D[pointsI][1] =
            pointsIn3DTranslated[pointsI] & orthoBase[1];

        // w component = error above projection plane
        dist[pointsI] = pointsIn3DTranslated[pointsI] & orthoBase[2];
    }

    distanceProjection = dist;

    return pointsIn2D;
}


// Projection of 2D points onto a 1D normalized direction vector
template<class MasterPatch, class SlavePatch>
scalarField newGGIInterpolation<MasterPatch, SlavePatch>::projectPoints2Dto1D
(
    const vector2D&      normalizedProjectionDir,
    const point2D&       normalizedProjectionDirOffset,
    const List<point2D>& lPoints2D
    ) const
{
    scalarField pointsIn1D(lPoints2D.size()); // Return values.

    forAll (lPoints2D, ptsI)
    {
        pointsIn1D[ptsI] =
            normalizedProjectionDir & (lPoints2D[ptsI]
          - normalizedProjectionDirOffset);
    }

    return pointsIn1D;
}


// Polygon overlap or polygon collision detection using the Separating
// Axes Theorem.  We expect this algorithm to possibly call himself a
// second time: code needs to be re-entrant!!!
//
// First time it is called, we test polygon2 agains each edges of
// polygon1 If we still detect an overlap, then the algorith will call
// itself to test polygon1 against polygon2 in order to find a possibe
// separating Axes.  The flag firstCall will be used for detecting if
// we are on the first or second call, so we can exit properly after
// the second call.
template<class MasterPatch, class SlavePatch>
bool newGGIInterpolation<MasterPatch, SlavePatch>::detect2dPolygonsOverlap
(
    const List<point2D>& poly1,
    const List<point2D>& poly2,
    const scalar& tolFactor,
    const bool firstCall
) const
{
    bool isOverlapping = true;

    for
    (
        label istart1 = 0, iend1 = poly1.size() - 1;
        istart1 < poly1.size();
        iend1 = istart1, istart1++
    )
    {
        // Reference edge from polygon1
        point2D curEdge1 = poly1[istart1] - poly1[iend1];

        point2D origEdge1 = poly1[iend1];

        // normalPerpDirection1 & curEdge1 == 0
        vector2D normalPerpDirection1(curEdge1[1], - curEdge1[0]);

        // Normalize
        normalPerpDirection1 /= mag(normalPerpDirection1) + VSMALL;

        // Project poly1 onto normalPerpDirection1
        scalarField poly1PointsProjection = projectPoints2Dto1D
        (
            normalPerpDirection1,
            origEdge1,
            poly1
        );

        // Project poly2 onto normalPerpDirection1
        scalarField poly2PointsProjection = projectPoints2Dto1D
        (
            normalPerpDirection1,
            origEdge1,
            poly2
        );

        // Grab range of the projected polygons
        scalar p1_min = Foam::min(poly1PointsProjection);
        scalar p1_max = Foam::max(poly1PointsProjection);
        scalar p2_min = Foam::min(poly2PointsProjection);
        scalar p2_max = Foam::max(poly2PointsProjection);

        // Here are all the possible situations for the range overlap
        // detection, and a simple test that discriminates them all
        //
        //
        //    P1 -------------                       p1_min + epsilon < p2_max  &&  p2_min + epsilon < p1_max  == true
        //    P2 -------------
        //
        //    P2 -------------                       p1_min + epsilon < p2_max  &&  p2_min + epsilon < p1_max  == true
        //    P1 -------------
        //
        //
        //    P1  -------------                      p1_min + epsilon < p2_max  &&  p2_min + epsilon < p1_max  == true
        //         P2  ----
        //
        //         P1  ----
        //    P2  -------------                      p1_min + epsilon < p2_max  &&  p2_min + epsilon < p1_max  == true
        //
        //
        //    P1  -------------                      p1_min + epsilon < p2_max  &&  p2_min + epsilon < p1_max  == true
        //              P2  ---------
        //
        //    P2  -------------                      p1_min + epsilon < p2_max  &&  p2_min + epsilon < p1_max  == true
        //              P1  ---------
        //
        //
        //                        P1 -------------   p1_min + epsilon < p2_max  &&  p2_min + epsilon < p1_max  == false
        //    P2  -------------
        //
        //
        //                        P2 -------------   p1_min + epsilon < p2_max  &&  p2_min + epsilon < p1_max  == false
        //    P1  -------------
        //
        //
        //
        //    P1  -------------                      p1_min + epsilon < p2_max  &&  p2_min + epsilon < p1_max  == false
        //                     ------------- P2
        //
        //
        //    P2  -------------                      p1_min + epsilon < p2_max  &&  p2_min + epsilon < p1_max  == false
        //                     ------------- P1
        //
        //
        //
        // What value should we give to epsilon???
        //
        // The intervals [p1_min, p1_max] and [p2_min, p2_max] are
        // basically the dimension of one side of the BB for a given
        // polygon, for a given "orientation" of the polygon.
        //
        // This means that epsilon^2 is roughly the size of the minimal
        // surface area intersecting the 2 polygons that one accept to
        // discard.
        //
        // So, if for instance we fix epsilon = 10e-3 * min
        // (range_of_polygon1, range_of_polygon2), this means that we
        // accept to discard an intersecting area roughly 10e-6 times the
        // surface of the smallest polygon, so 1 PPM.
        //
        // Which means also that our GGI weighting factors will never be
        //  smaller than roughly 10e-6, because this is the fraction of
        // the intersection surface area we choose to discard.

        scalar _epsilon = tolFactor*
            Foam::min
            (
                mag(p1_max - p1_min),
                mag(p2_max - p2_min)
            );

        // Evaluate the presence or not of an overlapping
        if
        (
            !((p1_min + _epsilon < p2_max)
          && (p2_min + _epsilon < p1_max))
        )
        {
            isOverlapping = false; // We have found a separating axis!
            break;
        }
    }

    if (isOverlapping && firstCall)
    {
        // We have not found any separating axes by exploring from
        // poly1, let's switch by exploring from poly2 instead
        isOverlapping =
            detect2dPolygonsOverlap(poly2, poly1, tolFactor, false);
    }

    return isOverlapping;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
