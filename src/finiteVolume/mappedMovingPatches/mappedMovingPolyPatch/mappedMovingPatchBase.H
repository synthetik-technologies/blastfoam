/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

Class
    Foam::mappedMovingPatchBase

Description
    Determines a mapping between patch face centres and mesh cell or face
    centres and processors they're on.

    If constructed from dictionary:
    \verbatim
        // Region to sample (default is region0)
        sampleRegion region0;

        // What to sample:
        // - nearestCell         : sample cell containing point
        // - nearestOnlyCell     : nearest sample cell (even if not containing
        //                         point)
        // - nearestPatchFace    : nearest face on selected patch
        // - nearestPatchFaceAMI : nearest face on selected patch
                                   - patches need not conform
                                   - uses AMI interpolation
        // - nearestFace         : nearest boundary face on any patch
        // - nearestPatchPoint   : nearest patch point (for coupled points
        //                         this might be any of the points so you have
        //                         to guarantee the point data is synchronised
        //                         beforehand)
        sampleMode nearestCell;

        // If sampleMode is nearestPatchFace : patch to find faces of
        samplePatch movingWall;

        // If sampleMode is nearestPatchFace : specify patchgroup to find
        // samplePatch and sampleRegion (if not provided)
        coupleGroup baffleGroup;

        // How to supply offset (w.r.t. my patch face centres):
        // - uniform : single offset vector
        // - nonuniform : per-face offset vector
        // - normal : using supplied distance and face normal
        offsetMode uniform;

        // According to offsetMode (see above) supply one of
        // offset, offsets or distance
        offset  (1 0 0);
    \endverbatim

    Note: if offsetMode is \c normal it uses outwards pointing normals. So
    supply a negative distance if sampling inside the domain.


Note
    Storage is not optimal. It temporary collects all (patch)face centres
    on all processors to keep the addressing calculation simple.

SourceFiles
    mappedMovingPatchBase.C

\*---------------------------------------------------------------------------*/

#ifndef mappedMovingPatchBase_H
#define mappedMovingPatchBase_H

#include "pointField.H"
#include "Tuple2.H"
#include "pointIndexHit.H"
#include "AMIInterpolation.H"
#include "coupleGroupIdentifier.H"
#include "fvMesh.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyPatch;
class polyMesh;
class mapDistribute;

/*---------------------------------------------------------------------------*\
                       Class mappedMovingPatchBase Declaration
\*---------------------------------------------------------------------------*/

class mappedMovingPatchBase
{

public:

    //- Helper class for finding nearest
    //  Nearest:
    //  - point+local index
    //  - sqr(distance)
    //  - processor
    typedef Tuple2<pointIndexHit, Tuple2<scalar, label>> nearInfo;

    class nearestEqOp
    {

    public:

        void operator()(nearInfo& x, const nearInfo& y) const
        {
            if (y.first().hit())
            {
                if (!x.first().hit())
                {
                    x = y;
                }
                else if (y.second().first() < x.second().first())
                {
                    x = y;
                }
            }
        }
    };

    class maxProcEqOp
    {

    public:

        void operator()(nearInfo& x, const nearInfo& y) const
        {
            if (y.first().hit())
            {
                if (!x.first().hit())
                {
                    x = y;
                }
                else if (y.second().second() > x.second().second())
                {
                    x = y;
                }
            }
        }
    };


protected:

    // Protected data

        //- Patch to sample
        const polyPatch& patch_;

        //- Region to sample
        mutable word sampleRegion_;


        //- Patch (if in sampleMode NEARESTPATCH*)
        mutable word samplePatch_;

        //- PatchGroup (if in sampleMode NEARESTPATCH*)
        const coupleGroupIdentifier coupleGroup_;

        //- Pointer to saved offset vector
        const volVectorField* DPtr_;


        // Derived information

            //- Communication schedule:
            //
            //    - Cells/faces to sample per processor
            //    - Patch faces to receive per processor
            //    - schedule
            mutable autoPtr<mapDistribute> mapPtr_;


    // Protected Member Functions

        //- Get the points from face-centre-decomposition face centres
        //  and project them onto the face-diagonal-decomposition triangles.
        tmp<pointField> facePoints(const polyPatch&) const;

        //- Collect single list of samples and originating processor+face.
        void collectSamples
        (
            const pointField& facePoints,
            pointField&,
            labelList& patchFaceProcs,
            labelList& patchFaces,
            pointField& patchFc
        ) const;

        //- Find cells/faces containing samples
        void findSamples
        (
            const pointField&,
            labelList& sampleProcs,     // processor containing sample
            labelList& sampleIndices,   // local index of cell/face
            pointField& sampleLocations // actual representative location
        ) const;

        //- Calculate mapping
        void calcMapping() const;


public:

    //- Runtime type information
    TypeName("mappedMovingPatchBase");


    // Constructors

        //- Construct from patch
        mappedMovingPatchBase(const polyPatch&);


        //- Construct from offsetMode=uniform
        mappedMovingPatchBase
        (
            const polyPatch& pp,
            const word& sampleRegion,
            const word& samplePatch
        );

        //- Construct from dictionary
        mappedMovingPatchBase(const polyPatch&, const dictionary&);

        //- Construct as copy, resetting patch
        mappedMovingPatchBase(const polyPatch&, const mappedMovingPatchBase&);

        //- Construct as copy, resetting patch, map original data
        mappedMovingPatchBase
        (
            const polyPatch&,
            const mappedMovingPatchBase&,
            const labelUList& mapAddressing
        );


    //- Destructor
    virtual ~mappedMovingPatchBase();


    // Member Functions

        void clearOut();

        // Access

            //- Region to sample
            inline const word& sampleRegion() const;

            //- Patch (only if NEARESTPATCHFACE)
            inline const word& samplePatch() const;

            //- PatchGroup (only if NEARESTPATCHFACE)
            inline const word& coupleGroup() const;

            //- Return size of mappedMoving mesh/patch/boundary
            inline label sampleSize() const;

            //- Set pointer to offset field
            inline void setOffsets(const volVectorField& D);

            //- Return reference to the parallel distribution map
            inline const mapDistribute& map() const;

            //- Get the region mesh
            const polyMesh& sampleMesh() const;

            //- Get the patch on the region
            const polyPatch& samplePolyPatch() const;


        // Helpers

            //- Get the sample points
            tmp<pointField> samplePoints() const;

            //- Get a point on the face given a face decomposition method:
            //  face-centre-tet : face centre. Returns index of face.
            //  face-planes     : face centre. Returns index of face.
            //  face-diagonal   : intersection of ray from cellcentre to
            //                    facecentre with any of the triangles.
            //                    Returns index (0..size-2) of triangle.
            static pointIndexHit facePoint
            (
                const polyMesh&,
                const label facei,
                const polyMesh::cellDecomposition
            );


        // Distribute

            //- Wrapper around map/interpolate data distribution
            template<class Type>
            void distribute(List<Type>& lst) const;

            //- Wrapper around map/interpolate data distribution with operation
            template<class Type, class CombineOp>
            void distribute(List<Type>& lst, const CombineOp& cop) const;

            //- Wrapper around map/interpolate data distribution
            template<class Type>
            void reverseDistribute(List<Type>& lst) const;

            //- Wrapper around map/interpolate data distribution with operation
            template<class Type, class CombineOp>
            void reverseDistribute(List<Type>& lst, const CombineOp& cop) const;


        // I/O

            //- Write as a dictionary
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mappedMovingPatchBaseI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mappedMovingPatchBaseTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
