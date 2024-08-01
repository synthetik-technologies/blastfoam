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
#include "primitivePatchInterpolation.H"

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
        const polyPatch& polyPatch_;

        //- Optional displacement field
        word displacementField_;

        //- Use inverse displacement
        bool inverseDisplacement_;

        // Demand-driven private data

            //- Displaced patch points
            mutable autoPtr<pointField> displacedPointsPtr_;

            //- Displaced/local patch
            mutable autoPtr<primitivePatch> physicalPatchPtr_;

            //- Old displaced patch points
            mutable autoPtr<pointField> displacedPoints0Ptr_;

            //- Old primitive patch
            mutable autoPtr<primitivePatch> physicalPatch0Ptr_;

            //- Point to face interpolator
            mutable autoPtr<primitivePatchInterpolation>
                pointToFaceInterpolatorPtr_;

            // Face to point weights
            mutable List<scalarField>* faceToPointWeightsPtr_;

            // Sum of face to point weights
            mutable scalarField* faceToPointSumWeightsPtr_;

    // Protected Member Functions

        //- Build global primitive patch
        virtual void calcPhysicalPatch() const;

        //- Make the point to face interpolator
        virtual void calcPointToFaceInterpolation() const;

        //- Make the face to point interpolator
        virtual void calcFaceToPointInterpolation() const;

        //- Return face to point interpolation weights
        const List<scalarField>& faceToPointWeights() const;

        //- Return sum of face to point interpolation weights
        const scalarField& faceToPointSumWeights() const;

        //- Return reference to the local interpolator
        const primitivePatchInterpolation& pointToFaceInterpolator() const;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        globalPolyPatch(const globalPolyPatch&);

        //- Disallow default bitwise assignment
        void operator=(const globalPolyPatch&);

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

        //- Construct from components
        globalPolyPatch
        (
            const polyPatch& patch,
            const word& displacementField = "none"
        );


    declareRunTimeSelectionTable
    (
        autoPtr,
        globalPolyPatch,
        patch,
        (
            const dictionary& dict,
            const polyPatch& patch
        ),
        (dict, patch)
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

        //- Return reference to patch
        const polyPatch& patch() const
        {
            return polyPatch_;
        }

        //- Set displacement field name
        void setDisplacementField(const word& Dname)
        {
            displacementField_ = Dname;
        }

        //- Set inverse displacement flag
        void setInverseDisplacement(const bool inv)
        {
            inverseDisplacement_ = inv;
        }

        //- Return if the global patch is currently set
        bool valid() const
        {
            return physicalPatchPtr_.valid();
        }

        //- Return the physical patch
        const primitivePatch& physicalPatch() const;

        //- Return the old physical patch
        const primitivePatch& physicalPatch0() const;


        // Interpolation

            //- Interpolate patch faces to points
            //  local face -> global face -> global point -> local point
            template<class Type>
            tmp<Field<Type>> faceToPoint(const Field<Type>& fFld) const;

            template<class Type>
            tmp<Field<Type>> faceToPoint(const tmp<Field<Type>>& tfFld) const;

            //- Interpolate patch points to faces
            //  points do not need to be made global
            template<class Type>
            tmp<Field<Type>> pointToFace(const Field<Type>& pFld) const;

            template<class Type>
            tmp<Field<Type>> pointToFace(const tmp<Field<Type>>& tpFld) const;


        //- Force updat
        virtual void update();

        //- Correct patch after moving points
        virtual void movePoints(const bool clear);

        //- Update for changes in topology
        virtual void updateMesh();

        //- Clear addressing
        virtual void clearOut() const;

        //- Write the global patch
        virtual bool write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "globalPolyPatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
