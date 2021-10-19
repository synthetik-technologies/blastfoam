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

Class
    newGGIInterpolation

Description
    Mass-conservative face interpolation of face data between two
    primitivePatches

    Note on parallelisation
    GGI search and cutting is an expensive operation, typically used in
    patch-to-patch interpolation.  The current paradigm keeps primitive patch
    data as a global zone on all processors for ease of manipulation
    and to avoid global numbering.

    newGGIInterpolation uses a globalData flag to indicate that identical
    patch data is available everywhere.  In such cases, ALL processors
    (not just the ones which hold a piece of the active GGI surface)
    can be used for the cutting, with the resulting addressing and weights
    data assembled by global communication at the end.

    If the GGI patch data is not identical on all processors, set
    globalData to false.  HJ, 27/Apr/2016

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

Contributor:
    Martin Beaudoin, Hydro-Quebec, (2008)

SourceFiles
    newGGIInterpolation.C
    newGGIInterpolate.C
    newGGIInterpolationWeights.C

\*---------------------------------------------------------------------------*/

#ifndef newGGIInterpolation_H
#define newGGIInterpolation_H

#include "className.H"
#include "labelList.H"
#include "scalarField.H"
#include "pointField.H"
#include "FieldFields.H"
#include "faceList.H"
#include "intersection.H"
#include "point2D.H"
#include "NamedEnum.H"
#include "tolerancesSwitch.H"
#include "optimisationSwitch.H"
#include "labelPair.H"
#include "boundBox.H"
#include "triPointRef.H"
#include "Map.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class newGGIInterpolationName Declaration
\*---------------------------------------------------------------------------*/

class newGGIInterpolationName
{
public:

    // Public enumerations

        //- Define a return type for isVertexInsidePolygon.
        //  The names are rather obvious
        // on vertex, or on edge or inside == inside
        enum insideOutside
        {
            ALL_OUTSIDE,
            ALL_INSIDE,
            PARTIALLY_OVERLAPPING
        };

        //- Quick reject overlap search algorithm
        enum quickReject
        {
            THREE_D_DISTANCE,
            AABB,
            BB_OCTREE
        };


    // Static data

        ClassName("newGGIInterpolation");

        //- Quick reject names
        static const NamedEnum<quickReject, 3> quickRejectNames_;


    // Constructors

        //- Construct null
        newGGIInterpolationName()
        {}

};


/*---------------------------------------------------------------------------*\
                      Class newGGIInterpolation Declaration
\*---------------------------------------------------------------------------*/

template<class MasterPatch, class SlavePatch>
class newGGIInterpolation
:
    public newGGIInterpolationName
{
    // Private data types

        // Definition for a 3D orthoNormalBasis type
        typedef VectorSpace<Vector<vector>, vector, 3> orthoNormalBasis;


    // Private data

        //- Reference to the master patch
        const MasterPatch& masterPatch_;

        //- Reference to the slave patch
        const SlavePatch& slavePatch_;

        //- Slave-to-master transformation tensor.  Transforms slave data to
        //  master plane.  Size equals number of slave faces; zero length
        //  indicates no transform.  Size 1 indicates constant transform
        //  HJ, 6/Jan/2009
        tensorField forwardT_;

        //- Master-to-slave transformation tensor.  Transforms slave data to
        //  master plane.  Size equals number of master faces; zero length
        //  indicates no transform.  Size 1 indicates constant transform
        //  HJ, 6/Jan/2009
        tensorField reverseT_;

        //- Slave-to-master separation vector. Translation of slave data to
        //  master plane.  Size equals number of slave faces; zero length
        //  indicates no translation.  MB, 28/Jan/2009
        vectorField forwardSep_;

        //- Global data: are master and slave present on all processors?
        //  If true, cutting will be performed in parallel
        bool globalData_;

        //- Master non-overlap face tolerance factor
        const scalar masterNonOverlapFaceTol_;

        //- Slave non-overlap face tolerance factor
        const scalar slaveNonOverlapFaceTol_;

        //- Rescale the GGI weighting factors or not
        const bool rescaleGGIWeightingFactors_;

        //- Quick reject search algorithm
        quickReject reject_;

        //- Use the previous candidate master neighbors in the quick reject
        //  search
        Switch usePrevCandidateMasterNeighbors_;

        //- Previous candidate master neighbors
        mutable labelListList prevCandidateMasterNeighbors_;

        //- Optional: region of interest where weights and intersections are
        //  calculated; outside this box, the intersections are not caculated.
        //  This is motivated by contact simulations where in general we might
        //  know the approximate contacting region, so there is no need to
        //  spend time checking faces outside of this
        const boundBox regionOfInterest_;


    // Demand-driven data

        // Master-to-slave interpolation

            //- Master to slave addressing
            mutable labelListList* masterAddrPtr_;

            //- Master to slave  weights
            mutable scalarListList* masterWeightsPtr_;

            //- Slave face into which each point of master patch is projected
            mutable List<labelPair>* masterPointAddressingPtr_;

            //- Weighting factors
            mutable FieldField<Field, scalar>* masterPointWeightsPtr_;

            //- Distance to intersection for master patch points
            mutable scalarField* masterPointDistancePtr_;

            //- Distance to intersection vectors for master patch points
            mutable vectorField* masterPointDistanceVectorsPtr_;

            //- Master edge loops map
            mutable List< Map<label> > masterEdgeLoopsMap_;

        // Slave-to-master interpolation

            //- Slave to master face addressing
            mutable labelListList* slaveAddrPtr_;

            //- Slave to master face weights
            mutable scalarListList* slaveWeightsPtr_;

            //- Master face into which each point of slave patch is projected
            mutable List<labelPair>* slavePointAddressingPtr_;

            //- Weighting factors
            mutable FieldField<Field, scalar>* slavePointWeightsPtr_;

            //- Distance to intersection for slave patch points
            mutable scalarField* slavePointDistancePtr_;

            //- Distance to intersection vector for slave patch points
            mutable vectorField* slavePointDistanceVectorsPtr_;

            //- Slave edge loops map
            mutable List< Map<label> > slaveEdgeLoopsMap_;

            //- Use the new method for calculating the point distances
            Switch useNewPointDistanceMethod_;

            //- Project points outside of contact back to the edge of the target
            //  patch
            Switch projectPointsToPatchBoundary_;

            //- Check and fix point distance orientations
            Switch checkPointDistanceOrientations_;

            //- List of uncovered master faces
            mutable labelList* uncoveredMasterAddrPtr_;

            //- List of uncovered slave faces
            mutable labelList* uncoveredSlaveAddrPtr_;

  // Private static data

        //- Facet area error tolerance
        static const debug::tolerancesSwitch areaErrorTol_;

        //- Facet normal featureCos criteria
        static const debug::tolerancesSwitch featureCosTol_;

        //- Facet bound box extension factor
        static const debug::tolerancesSwitch faceBoundBoxExtendSpanFraction_;

       //- The next 3 attributes are parameters controlling the creation
       //  of an octree search engine for the GGI facets neighbourhood
       //  determination.
       //
       //  Tests using GGI patches of up to ~100K facets have shown that
       //  the following values gave the best compromise between the
       //  "octree creation time" vs "octree search time":
       //
       //     octreeSearchMinNLevel_     : 3
       //     octreeSearchMaxLeafRatio_  : 3
       //     octreeSearchMaxShapeRatio_ : 1
       //
       //  For GGI patches larger than ~100K facets, your mileage may vary.
       //  So these 3 control parameters are adjustable using the following
       //  global optimization switches:
       //
       //     GGIOctreeSearchMinNLevel
       //     GGIOctreeSearchMaxLeafRatio
       //     GGIOctreeSearchMaxShapeRatio
       //

       //- Octree search: minNlevel parameter for octree constructor
       static const debug::optimisationSwitch octreeSearchMinNLevel_;

       //- Octree search: maxLeafRatio parameter for octree constructor
       static const debug::optimisationSwitch octreeSearchMaxLeafRatio_;

       //- Octree search: maxShapeRatio parameter for octree constructor
       static const debug::optimisationSwitch octreeSearchMaxShapeRatio_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        newGGIInterpolation(const newGGIInterpolation&);

        //- Disallow default bitwise assignment
        void operator=(const newGGIInterpolation&);

        //- Calculate point weights
        void calcMasterPointAddressing() const;

        //- Calculate point weights
        void calcMasterPointWeights() const;

        //- Calculate point weights
        void calcSlavePointAddressing() const;

        //- Calculate point weights
        void calcSlavePointWeights() const;

        //- Calculate edge loops map
        void calcMasterEdgeLoopsMap() const;

        //- Calculate edge loops map
        void calcSlaveEdgeLoopsMap() const;

        //- Return edge loops map
        const List< Map<label> >& masterEdgeLoopsMap() const;

        //- Return edge loops map
        const List< Map<label> >& slaveEdgeLoopsMap() const;

    // Helper functions for parallel search

        //- Is parallel cutting active?
        bool globalData() const
        {
            return
                Pstream::parRun()
             && globalData_
             && masterPatch_.size() > 10*Pstream::nProcs();
        }

        //- Parallel master slice start
        label parMasterStart() const;

        //- Parallel master slice end
        label parMasterEnd() const;

        //- Parallel master slice size
        label parMasterSize() const;


    // Helper functions for demand-driven data

        //- Evaluate faces neighborhood based of sphere defined by faces BB
        void findNeighbours3D(labelListList& result) const;

        //- Evaluate faces neighborhood based of faces Axis Aligned BB
        void findNeighboursAABB(labelListList& result) const;

        //- Evaluate faces neighborhood based of faces BB and an octree
        //  search engine
        void findNeighboursBBOctree(labelListList& result) const;

        //- Update the neighbour faces using the previous neighbours and the
        //  Axis Aligned BB method
        void updateNeighboursAABB(labelListList& result) const;

        //- Projects a list of points onto a plane located at
        //  planeOrig, oriented along planeNormal
        tmp<pointField> projectPointsOnPlane
        (
            const pointField& lpoints,
            const vector& planeOrig,
            const vector& planeDirection,
            scalarField& distanceProjection
        ) const;

        //- Compute an orthonormal basis. Useful for 2D projection
        orthoNormalBasis computeOrthonormalBasis
        (
            const vector& normalVectorCentre,
            const vector& normalVector,
            const pointField& pointsOnPlane
        ) const;

        //- Compute the 2D projection of a 3D list of point into the
        //  orthonormal basis orthoBase
        List<point2D> projectPoints3Dto2D
        (
            const orthoNormalBasis& orthoBase,
            const vector& orthoBaseOffset,
            const pointField& pointsIn3D,
            scalarField& distanceProjection
        ) const;

        //- Compute the 1D projection of a 2D list of points onto the
        //  projectionDir direction
        scalarField projectPoints2Dto1D
        (
            const vector2D& normalizedProjectionDir,
            const point2D& normalizedProjectionDirOffset,
            const List<point2D>& lPoints2D
        ) const;

        //- Detect the intersection of 2 2D polygons using the method of
        //  Separating Axes Theorem
        bool detect2dPolygonsOverlap
        (
            const List<point2D>& poly1,
            const List<point2D>& poly2,
            const scalar& tolFactor,
            const bool firstCall = true // Optional parameter: first call
        ) const;

        //- Compute the intersection area between 2 2D polygons
        scalar polygonIntersection
        (
            const List<point2D>& poly1,
            const List<point2D>& poly2
        ) const;

        //- Detect is the points from a subject polygon are inside or
        //  outside a clipping polygon
        insideOutside isVertexInsidePolygon
        (
            const List<point2D>& clippingPolygon,
            const List<point2D>& subjectPolygon,
            List<bool>& subjectVertexInside
        ) const;

        //- Compute the intersection between two 2D polygons using the
        //  Sutherland Hodgman algorithm
        List<point2D> clipPolygon2DSutherlandHodgman
        (
            const List<point2D>& clippingPolygon,
            const List<point2D>& subjectPolygon
        ) const;

        //- Compute the intersection between two 2D polygons using the
        //  Greiner Hormann algorithm
        List<point2D> clipPolygon2DGreinerHormann
        (
            const List<point2D>& clippingPolygon,
            const List<point2D>& subjectPolygon,
            const List<bool>& subjectVertexInside
        ) const;

        //- Compute the area of a 2D polygon
        scalar area2D
        (
            const List<point2D>& polygon
        ) const;


        //- Calculate addressing and weights
        void calcAddressing() const;

        //- Rescale GGI weighting factors
        void rescaleWeightingFactors() const;

        //- Find non-overlapping faces
        tmp<labelField> findNonOverlappingFaces
        (
            const scalarListList& patchWeights,
            const scalar& nonOverlapFaceTol
        ) const;

        //- Set the source point distance, vcetor distance and addressing
        void setSourceToTargetPointAddressing
        (
            List<labelPair>& sourcePointAddr, // source point addressing
            scalarField& sourcePointDist,     // source point distances
            vectorField& sourcePointDistVecs, // source point distance vectors
            const MasterPatch& sourcePatch,   // source patch
            const SlavePatch& targetPatch,    // target patch
            const labelListList& sourceAddr,   // source patch face addressing
            const List< Map<label> >& targetEdgeLoopsMap // target edge loop map
        ) const;

        //- Set the source point distance and addressing using the previous
        //  method
        void setSourceToTargetPointAddressingPrev
        (
            List<labelPair>& sourcePointAddr, // source point addressing
            scalarField& sourcePointDist,     // source point distances
            const MasterPatch& sourcePatch,   // source patch
            const SlavePatch& targetPatch,    // target patch
            const labelListList& sourceAddr   // source patch face addressing
        ) const;

        //- Return list of possibly in contact shadow faces for the given faces
        tmp<labelField> possibleTargetFaces
        (
            const labelList& sourceFaces,
            const labelListList& sourceFaceAddr
        ) const;

        //- Calculates the list of triangular sub-faces to which the current
        //  point can be projected.
        //  Each labelPair entry returns the master face index (first) and the
        //  face local point index indicating the sub-triangle
        //  The projection directions are the non-unit vectors from P to the
        //  triangular sub face
        void projectPointToTriangularSubFaces
        (
            List<labelPair>& triSubFacesIndex,
            List<vector>& triSubFacesProjDirections,
            PtrList< triangle<point, point> >& triSubFaces,
            const point& P,
            const labelList& possibleTargetFaces,
            const faceList& targetFaces,
            const pointField& taretPoints
        ) const;

        //- Find the closest vertex to P in the given set of faces
        label findClosestVertex
        (
            const point& P,
            const labelList& possibleTargetFaces,
            const faceList& targetFaces,
            const pointField& targetPoints
        ) const;

        //- Do the triangles share an edge
        void triangleSharedPoints
        (
            List<point>& sharedPoints,
            const triangle<point, point>& triA,
            const triangle<point, point>& triB
        ) const;

        //- Correct point projections on the boundary of the patch
        void correctBoundaryPointProjection
        (
            const MasterPatch& targetPatch,
            const label closestTargetPointID,
            const labelListList& targetEdgeLoops,
            const List< Map<label> >& targetEdgeLoopsMap,
            const vector& MP,
            scalar& sourcePointDist,
            vector& sourcePointDistVecs
        ) const;

        //- Find the boundary point bi-normal
        //  This is the vector pointing away from the patch internal field
        void boundaryPointBiNormal
        (
            vector& boundaryPointBiNormal,
            const MasterPatch& targetPatch,
            const label closestTargetPointID,
            const labelList& targetBoundaryPoints,
            const Map<label>& targetBoundaryPointsMap
        ) const;

        //- Clear all geometry and addressing
        void clearOut();


    // Interpolation data access

        //- Interpolate given source and target, addressing and weights
        template<class Type>
        static void interpolate
        (
            const Field<Type>& ff,
            Field<Type>& result,
            const labelListList& addr,
            const scalarListList& weights
        );

        //- Interpolate given source and target, addressing and weights
        //  for masked faces only
        template<class Type>
        static void maskedInterpolate
        (
            const Field<Type>& ff,
            Field<Type>& result,
            const labelList& mask,
            const labelListList& addr,
            const scalarListList& weights
        );

        //- Bridge uncovered faces given addressing
        template<class Type>
        static void bridge
        (
            const Field<Type>& bridgeField,
            Field<Type>& ff,
            const labelList& addr
        );

        //- Bridge uncovered faces given addressing
        //  for masked faces only
        template<class Type>
        static void maskedBridge
        (
            const Field<Type>& bridgeField,
            Field<Type>& ff,
            const labelList& mask,
            const labelList& uncoveredFaces
        );

        //- Is a transform required?
        inline bool doTransform() const
        {
            return forwardT_.size() > 0;
        }

        //- Is a separation required?
        inline bool doSeparation() const
        {
            return forwardSep_.size() > 0;
        }


public:

    // Constructors

        //- Construct from components
        newGGIInterpolation
        (
            const MasterPatch& masterPatch,
            const SlavePatch&  slavePatch,
            const tensorField& forwardT,
            const tensorField& reverseT,
            const vectorField& forwardSep,
            const bool globalData,
            const scalar masterFaceNonOverlapFaceTol = 0,
            const scalar slaveFaceNonOverlapFaceTol = 0,
            const bool rescaleGGIWeightingFactors = true,
            const quickReject reject = BB_OCTREE,
            const boundBox regionOfInterest =
                boundBox(vector::min, vector::max)
        );


    // Destructor
    ~newGGIInterpolation();


    // Member Functions

        // Access

            //- Reference to the master patch
            const MasterPatch& masterPatch() const
            {
                return masterPatch_;
            }

            //- Reference to the slave patch
            const SlavePatch& slavePatch() const
            {
                return slavePatch_;
            }

            //- Return reference to master addressing
            const labelListList& masterAddr() const;

            //- Return reference to master weights
            const scalarListList& masterWeights() const;

            //- Return reference to slave addressing
            const labelListList& slaveAddr() const;

            //- Return reference to slave weights
            const scalarListList& slaveWeights() const;

            //- Return reference to the master list of non-overlap faces
            const labelList& uncoveredMasterFaces() const;

            //- Return reference to the slave list of non-overlap faces
            const labelList& uncoveredSlaveFaces() const;

            //- Return reference to point addressing
            const List<labelPair>& masterPointAddr() const;

            //- Return reference to point weights
            const FieldField<Field, scalar>& masterPointWeights() const;

            //- Return distance to intersection for patch points
            const scalarField& masterPointDistanceToIntersection() const;

            //- Return distance to intersection vectors for patch points
            const vectorField& masterPointDistanceVectorsToIntersection() const;

            //- Return reference to point addressing
            const List<labelPair>& slavePointAddr() const;

            //- Return reference to point weights
            const FieldField<Field, scalar>& slavePointWeights() const;

            //- Return distance to intersection for patch points
            const scalarField& slavePointDistanceToIntersection() const;

            //- Return distance to intersection vectors for patch points
            const vectorField& slavePointDistanceVectorsToIntersection() const;

            //- Non-const reference to the usePrevCandidateMasterNeighbours
            //  switch
            Switch& usePrevCandidateMasterNeighbors()
            {
                return usePrevCandidateMasterNeighbors_;
            }

            //- Clear the prevCandidateMasterNeighbours list
            void clearPrevCandidateMasterNeighbors()
            {
                prevCandidateMasterNeighbors_.clear();
            }

            //- Non-const reference to the point distance calculation switch
            Switch& useNewPointDistanceMethod()
            {
                return useNewPointDistanceMethod_;
            }

            //- Non-const reference to the projectPointsToPatchBoundary switch
            Switch& projectPointsToPatchBoundary()
            {
                return projectPointsToPatchBoundary_;
            }

            //- Non-const reference to the checkPointDistanceOrientations switch
            Switch& checkPointDistanceOrientations()
            {
                return checkPointDistanceOrientations_;
            }

    // Interpolation functions

        //- Interpolate from master to slave
        template<class Type>
        tmp<Field<Type> > masterToSlave(const Field<Type>& ff) const;

        template<class Type>
        tmp<Field<Type> > masterToSlave
        (
            const tmp<Field<Type> >& tff
        ) const;

        //- Interpolate from master to slave, only for marked master faces
        //  Addressing picks the faces from patch that are selected
        //  for collection into the result field
        template<class Type>
        void maskedMasterToSlave
        (
            const Field<Type>& ff,
            Field<Type>& result,
            const labelList& mask
        ) const;

        //- Interpolate from slave to master
        template<class Type>
        tmp<Field<Type> > slaveToMaster(const Field<Type>& ff) const;

        template<class Type>
        tmp<Field<Type> > slaveToMaster
        (
            const tmp<Field<Type> >& tff
        ) const;

        //- Interpolate from slave to master, only for marked slave faces
        //  Addressing picks the faces from patch that are selected
        //  for collection into the result field
        template<class Type>
        void maskedSlaveToMaster
        (
            const Field<Type>& ff,
            Field<Type>& result,
            const labelList& mask
        ) const;


        //- Bridge uncovered master patch field
        template<class Type>
        void bridgeMaster
        (
            const Field<Type>& bridgeField,
            Field<Type>& ff
        ) const;

        //- Bridge uncovered master patch field, only for marked master faces
        template<class Type>
        void maskedBridgeMaster
        (
            const Field<Type>& bridgeField,
            Field<Type>& ff,
            const labelList& mask
        ) const;

        //- Bridge uncovered slave patch field
        template<class Type>
        void bridgeSlave
        (
            const Field<Type>& bridgeField,
            Field<Type>& ff
        ) const;

        //- Bridge uncovered slave patch field, only for marked slave faces
        template<class Type>
        void maskedBridgeSlave
        (
            const Field<Type>& bridgeField,
            Field<Type>& ff,
            const labelList& mask
        ) const;

        //- Interpolate point field
        template<class Type>
        tmp<Field<Type> > masterToSlavePointInterpolate
        (
            const Field<Type>& pf
        ) const;

//         template<class Type>
//         tmp<Field<Type> > masterToSlavePointInterpolate
//         (
//             const tmp<Field<Type> >& tpf
//         ) const;

        template<class Type>
        tmp<Field<Type> > slaveToMasterPointInterpolate
        (
            const Field<Type>& pf
        ) const;

//         template<class Type>
//         tmp<Field<Type> > slaveToMasterPointInterpolate
//         (
//             const tmp<Field<Type> >& tpf
//         ) const;


    // Edit

        //- Correct weighting factors for moving mesh.
        bool movePoints
        (
            const tensorField& forwardT,
            const tensorField& reverseT,
            const vectorField& forwardSep
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "newGGIInterpolation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
