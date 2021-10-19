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
    Foam::globalPolyPatch

Description
    A mesh patch synced in parallel runs such that all faces are present
    on all processors.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved
    Modifications/additions by Philip Cardiff, UCD.  All rights reserved

SourceFiles
    globalPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef globalPolyPatch_H
#define globalPolyPatch_H

#include "typeInfo.H"
#include "dictionary.H"
#include "standAlonePatch.H"
#include "polyMesh.H"
#include "PrimitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
                      Class globalPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class globalPolyPatch
{
protected:
    // Protected data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Patch name
        word patchName_;

        //- Reference to patch
        const polyPatch& patch_;

        //- Optional displacement field
        word displacementField_;

        // Demand-driven private data

            //- Primitive patch made out of faces from parallel decomposition
            mutable autoPtr<standAlonePatch> globalPatchPtr_;

            //- Patch point to zone addressing
            mutable autoPtr<labelList> pointToGlobalAddrPtr_;

            //- Patch face to zone addressing
            mutable autoPtr<labelList> faceToGlobalAddrPtr_;

            //- Addressing into master processor global patch from current
            //  processor global patch
            mutable autoPtr<labelList> globalMasterToCurrentProcPointAddrPtr_;

            //- Patch interpolator
            mutable autoPtr<PrimitivePatchInterpolation<standAlonePatch>> interpPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        globalPolyPatch(const globalPolyPatch&);

        //- Disallow default bitwise assignment
        void operator=(const globalPolyPatch&);

        //- Build global primitive patch
        void calcGlobalPatch() const;

        // Make globalMasterToCurrentProcPointAddr
        void calcGlobalMasterToCurrentProcPointAddr() const;

        //- Make the interpolator
        void calcInterp() const;

        //- Check
        void check() const;


public:

    //- Runtime type information
    TypeName("globalPolyPatch");


    // Constructors

        //- Construct from components
        globalPolyPatch
        (
            const dictionary& dict,
            const polyPatch& patch
        );


    //- Return a new globalPolyPatch
    static autoPtr<globalPolyPatch> New
    (
        const dictionary& dict,
        const polyPatch& patch
    );


    //- Destructor
    virtual ~globalPolyPatch();


    // Member Functions

        //- Return patch name
        const word& patchName() const
        {
            return patchName_;
        }

        //- Return mesh reference
        const polyMesh& mesh() const;

        //- Map storing the local face index for every global face index.
        //  Used to find out the index of face in the zone from the known
        //   global face index.  If the face is not in the zone, returns -1
        label whichFace(const label globalFaceID) const;

        //- Return reference to patch
        const polyPatch& patch() const
        {
            return patch_;
        }

        //- Return reference to global patch
        const standAlonePatch& globalPatch() const;

        //- Return reference to the interpolator
        const PrimitivePatchInterpolation<standAlonePatch>&
        interpolator() const;

        //- Return patch point to global point addressing
        const labelList& pointToGlobalAddr() const;

        //- Return patch face to global face addressing
        const labelList& faceToGlobalAddr() const;

        // Addressing into master processor global patch from current processor
        // global patch
        const labelList& globalMasterToCurrentProcPointAddr() const;

        // Interpolation

            //- Interpolate patch to global patch point data
            template<class Type>
            tmp<Field<Type> > patchPointToGlobal
            (
                const Field<Type>& pField
            ) const;

            template<class Type>
            tmp<Field<Type> > patchPointToGlobal
            (
                const tmp<Field<Type>>& pField
            ) const;

            //- Filter global patch point data to patch
            template<class Type>
            tmp<Field<Type> > globalPointToPatch
            (
                const Field<Type>& gField
            ) const;

            template<class Type>
            tmp<Field<Type> > globalPointToPatch
            (
                const tmp<Field<Type>>& gField
            ) const;

            //- Interpolate patch to global patch face data
            template<class Type>
            tmp<Field<Type> > patchFaceToGlobal
            (
                const Field<Type>& pField
            ) const;

            template<class Type>
            tmp<Field<Type> > patchFaceToGlobal
            (
                const tmp<Field<Type>>& pField
            ) const;

            //- Filter global patch face data to patch
            template<class Type>
            tmp<Field<Type> > globalFaceToPatch
            (
                const Field<Type>& gField
            ) const;

            template<class Type>
            tmp<Field<Type> > globalFaceToPatch
            (
                const tmp<Field<Type>>& gField
            ) const;


        //- Correct patch after moving points
        virtual void movePoints();

        //- Update for changes in topology
        virtual void updateMesh();

        //- Clear addressing
        virtual void clearOut() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "globalPolyPatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
