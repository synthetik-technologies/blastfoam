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

Class
    ggiPolyPatch

Description
    Generalised grid interface (GGI) patch.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

Contributor
    Martin Beaudoin, Hydro-Quebec, (2008)

SourceFiles
    ggiPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef ggiPolyPatch_H
#define ggiPolyPatch_H

#include "seperatedPolyPatch.H"
#include "ggiInterpolation.H"
#include "transform.H"
#include "word.H"
#include "faceZone.H"
#include "Switch.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ggiPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class ggiPolyPatch
:
    public seperatedPolyPatch
{
    // Private data

        //- Shadow patch name
        const word shadowName_;

        //- Interpolation zone name
        const word zoneName_;

        //- Use bridging to fix overlap error in interpolation
        Switch bridgeOverlap_;

        //- Quick reject algorithm
        ggiZoneInterpolation::quickReject reject_;


        // Demand-driven data

            //- Shadow patch index.  Delayed evaluation for construction
            mutable label shadowIndex_;

            //- Interpolation zone index.  Delayed evaluation for construction
            mutable label zoneIndex_;

            //- Patch-to-patch interpolation
            mutable ggiZoneInterpolation* patchToPatchPtr_;

            //- Patch-to-zone addressing
            mutable labelList* zoneAddressingPtr_;

            //- Remote zone addressing: data needed for interpolation
            mutable labelList* remoteZoneAddressingPtr_;

            //- Reconstructed patch neighbour cell centres
            mutable vectorField* reconFaceCellCentresPtr_;


        // Parallel communication optimisation, stored on master processor

            //- Is the patch localised on a single processor
            // (single processor in a parallel run)?
            //  Used for parallel optimisation
            mutable bool* localParallelPtr_;

            //- Communicator to use for parallel communication
            //  Note: only master side allocates the comm
            mutable label comm_;

            //- Message tag used for sending
            mutable int tag_;

            //- Map-distribute comms tool
            mutable mapDistribute* mapPtr_;


    // Private member functions

        //- Calculate patch-to-zone addressing
        virtual void calcZoneAddressing() const;

        //- Calculate remote patch-to-zone addressing
        virtual void calcRemoteZoneAddressing() const;

        //- Calculate interpolation
        virtual void calcPatchToPatch() const;


        // Geometry

            //- Calculate reconstructed cell centres
            void calcReconFaceCellCentres() const;

            //- Force calculation of transformation tensors
            virtual void calcTransforms() const;


        // Parallel communication optimisation, stored on master processor

            //- Calculate local parallel switch
            void calcLocalParallel() const;

            //- Calculate send and receive addressing
            void calcSendReceive() const;


        // Memory management

            //- Clear geometry
            void clearGeom();

            //- Clear out
            void clearOut();


protected:

    // Protected Member functions

        //- Is the GGI active? (zone and shadow present)
        bool active() const;

        //- Initialise the calculation of the patch addressing
//         virtual void initAddressing();

        //- Calculate the patch addressing
//         virtual void calcAddressing();

        //- Initialise the calculation of the patch geometry
        virtual void initCalcGeometry(PstreamBuffers&);

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&);

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&);

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&);


public:

    //- Runtime type information
    TypeName("ggi");


    // Constructors

        //- Construct from components
        ggiPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct from components
        ggiPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType,
            const word& shadowName,
            const word& zoneName,
            const bool bridgeOverlap,
            const ggiZoneInterpolation::quickReject
                reject = ggiZoneInterpolation::BB_OCTREE
        );

        //- Construct from dictionary
        ggiPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh&,
            const word& patchType
        );

        //- Construct as given the original patch and resetting the
        //  face list and boundary mesh information
        ggiPolyPatch
        (
            const ggiPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct as copy
        ggiPolyPatch(const ggiPolyPatch&) = delete;

        //- Construct as copy, resetting the boundary mesh
        ggiPolyPatch
        (
            const ggiPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new ggiPolyPatch
                (
                    refCast<const ggiPolyPatch>(*this),
                    bm,
                    index,
                    newSize,
                    newStart
                )
            );
        }


    //- Destructor
    virtual ~ggiPolyPatch();


    // Member functions

        // Basic info

            //- Return shadow patch name
            const word& shadowName() const
            {
                return shadowName_;
            }

            //- Return name of interpolation face zone
            const word& zoneName() const
            {
                return zoneName_;
            }

            //- Return shadow patch index
            label shadowIndex() const;

            //- Return zone patch index
            label zoneIndex() const;

            //- Return shadow patch
            const ggiPolyPatch& shadow() const;

            //- Return interpolation face zone
            const faceZone& zone() const;


        // Communications support

            //- Return communicator used for communication
            virtual label comm() const;

            //- Return message tag to use for communication
            virtual int tag() const;


        // Interpolation data

            //- Is this the master side?
            virtual bool owner() const
            {
                return index() < shadowIndex();
            }

            //- Return transformation between the coupled patches
            virtual const transformer& transform() const
            {
                return transformer::null;
            }

            //- Use bridging to fix overlap error in interpolation
            bool bridgeOverlap() const
            {
                return bridgeOverlap_;
            }

            //- Return patch-to-zone addressing
            const labelList& zoneAddressing() const;

            //- Return remote patch-to-zone addressing
            const labelList& remoteZoneAddressing() const;

            //- Is the patch localised on a single processor
            bool localParallel() const;

            //- Return mapDistribute
            const mapDistribute& map() const;

            //- Return reference to patch-to-patch interpolation
            //  Used only for addressing
            const ggiZoneInterpolation& patchToPatch() const;


        // Interpolation functions and bridging operations for fully uncovered
        // and partially covered faces

            //- Expand face field to full zone size, including reduction
            template<class Type>
            tmp<Field<Type> > fastExpand(const UList<Type>& pf) const;

            //- Interpolate face field: given field on a the shadow side,
            //  create an interpolated field on this side
            template<class Type>
            tmp<Field<Type> > interpolate(const Field<Type>& pf) const;

            template<class Type>
            tmp<Field<Type> > interpolate(const tmp<Field<Type> >& tpf) const;

            //- Set given field for fully uncovered faces
            template<class Type>
            void setUncoveredFaces
            (
                const Field<Type>& fieldToSet,
                Field<Type>& ff
            ) const;

            //- Set given field for partially covered faces
            template<class Type>
            void setPartialFaces
            (
                const Field<Type>& fieldToSet,
                Field<Type>& ff
            ) const;

            //- Scale field for partially covered faces
            template<class Type>
            void scalePartialFaces(Field<Type>& ff) const;

            //- Add given field to partially covered faces
            template<class Type>
            void addToPartialFaces
            (
                const Field<Type>& fieldToAdd,
                Field<Type>& ff
            ) const;


        // Geometric data

            //- Return reconstructed cell centres
            const vectorField& reconFaceCellCentres() const;


        // Patch ordering

            //- Initialize ordering for primitivePatch. Does not
            //  refer to *this (except for name() and type() etc.)
            virtual void initOrder(PstreamBuffers&, const primitivePatch&) const;

            //- Return new ordering for primitivePatch.
            //  Ordering is -faceMap: for every face
            //  index of the new face -rotation: for every new face the
            //  clockwise shift of the original face. Return false if nothing
            //  changes (faceMap is identity, rotation is 0), true otherwise.
            virtual bool order
            (
                PstreamBuffers&,
                const primitivePatch&,
                labelList& faceMap,
                labelList& rotation
            ) const;

            //- Synchronise communications of ordering for primitivePatch
            //  Used in cases when no topological change happens locally,
            //  but is happening on other processors
            virtual void syncOrder() const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ggiPolyPatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
