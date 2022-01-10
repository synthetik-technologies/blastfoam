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
    Foam::coupledGlobalPolyPatch

Description
    A mesh patch synced in parallel runs such that all faces are present
    on all processors.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved
    Modifications/additions by Philip Cardiff, UCD.  All rights reserved

SourceFiles
    coupledGlobalPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef coupledGlobalPolyPatch_H
#define coupledGlobalPolyPatch_H

#include "globalPolyPatch.H"
#include "standAlonePatch.H"
#include "polyMesh.H"
#include "patchToPatchMapping.H"
#include "UautoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
                      Class coupledGlobalPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class coupledGlobalPolyPatch
:
    public globalPolyPatch
{
private:
    // Private data

        //- Dictionary
        const dictionary& dict_;

        //- Region to sample
        word sampleRegion_;

        //- Patch to sample
        word samplePatch_;

        //- Unmapped faces
        mutable labelList unmappedFaces_;
        mutable labelList unmappedPoints_;


        // Demand-driven private data

            //- Patch interpolator
            mutable patchToPatchMapping*
                patchToPatchInterpPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        coupledGlobalPolyPatch(const coupledGlobalPolyPatch&);

        //- Disallow default bitwise assignment
        void operator=(const coupledGlobalPolyPatch&);

        //- Make the interpolator
        void calcPatchToPatchInterp() const;

        //- Set the interpolation pointer
        void setPatchToPatchInterp(patchToPatchMapping*) const;

        //- Clear the interpolator
        virtual void clearInterp(const bool top = true) const;


public:

    //- Runtime type information
    TypeName("coupledGlobalPolyPatch");


    // Constructors

        //- Construct from components
        coupledGlobalPolyPatch
        (
            const dictionary& dict,
            const polyPatch& patch
        );


    //- Destructor
    virtual ~coupledGlobalPolyPatch();


    // Member Functions

        //- Return patch name
        const word& sampleRegion() const
        {
            return sampleRegion_;
        }

        //- Return mesh reference
        const polyMesh& sampleMesh() const;

        //- Return reference to sampled globalPolyPatch
        const coupledGlobalPolyPatch& samplePatch() const;

        //- Return reference to the interpolator
        const patchToPatchMapping& patchToPatchInterpolator() const;

        //- Does the patch have unmapped faces
        bool hasUnmappedFaces() const
        {
            return unmappedFaces_.size();
        }

        //- Does the patch have unmapped points
        bool hasUnmappedPoints() const
        {
            return unmappedPoints_.size();
        }


        // Interpolation (from sample patch to this patch)
        // All return a field the size of this patch (face or points)

            //- Interpolate point field
            template<class Type>
            tmp<Field<Type>> pointInterpolate(const Field<Type>& pf) const;

            template<class Type>
            tmp<Field<Type>> pointInterpolate(const tmp<Field<Type>>& tpf) const;

            //- Interpolate face field
            template<class Type>
            tmp<Field<Type>> faceInterpolate(const Field<Type>& pf) const;

            template<class Type>
            tmp<Field<Type>> faceInterpolate(const tmp<Field<Type>>& tpf) const;

            //- Interpolate from faces to points
            template<class Type>
            tmp<Field<Type>> faceToPointInterpolate
            (
                const Field<Type>& ff
            ) const;

            template<class Type>
            tmp<Field<Type>> faceToPointInterpolate
            (
                const tmp<Field<Type>>& tff
            ) const;

            //- Interpolate from points to faces
            template<class Type>
            tmp<Field<Type>> pointToFaceInterpolate
            (
                const Field<Type>& pf
            ) const;

            template<class Type>
            tmp<Field<Type>> pointToFaceInterpolate
            (
                const tmp<Field<Type>>& tpf
            ) const;

            //- Set the unmapped values on field corresponding to faces
            template<class Type>
            void setUnmappedFace
            (
                Field<Type>& ff,
                const Type& unmapped
            ) const;

            template<class Type>
            void setUnmappedFace
            (
                Field<Type>& ff,
                const Field<Type>& unmapped
            ) const;

            template<class Type>
            void setUnmappedFace
            (
                Field<Type>& ff,
                const tmp<Field<Type>>& tunmapped
            ) const;

            //- Set the unmapped values on field corresponding to points
            template<class Type>
            void setUnmappedPoint
            (
                Field<Type>& pf,
                const Type& unmapped
            ) const;

            template<class Type>
            void setUnmappedPoint
            (
                Field<Type>& pf,
                const Field<Type>& unmapped
            ) const;

            template<class Type>
            void setUnmappedPoint
            (
                Field<Type>& pf,
                const tmp<Field<Type>>& tunmapped
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
#   include "coupledGlobalPolyPatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
