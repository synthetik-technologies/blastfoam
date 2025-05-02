/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::patchToPatch

Description
    Class to generate coupling geometry between two primitive patches

SourceFiles
    patchToPatch.C
    patchToPatchParallelOps.C
    patchToPatchTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef patchToPatchMapping_H
#define patchToPatchMapping_H

#include "runTimeSelectionTables.H"
#include "NamedEnum.H"
#include "dictionary.H"
#include "primitivePatch.H"
#include "treeBoundBox.H"
#include "standAlonePatch.H"
#include "remote.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class patchToPatchMapping Declaration
\*---------------------------------------------------------------------------*/

class patchToPatchMapping
{
protected:
    // Protected data

        //- Copy of settings dictionary
        const dictionary dict_;

        //- Is point mapping needed
        bool needPoints_;

        //- Stored patch pointers for determining mapping direction
        const primitivePatch& srcPatch_;
        const primitivePatch& tgtPatch_;

        //- Local data for parallel (non-singular processor)
        autoPtr<standAlonePatch> localTgtPatchPtr_;
        autoPtr<pointField> localTgtPoints0Ptr_;


        // Mapping data

            //- Are patches oriented in the same direction
            const bool reverse_;

            //- Single processor mapping
            label singleProcess_;


            // Face mapping

                //- Local source to target faces
                List<DynamicList<label>> localSrcFacesToTgt_;

                //- Local target to source faces
                List<DynamicList<label>> localTgtFacesToSrc_;

                //- For each source face, the coupled target weights
                List<DynamicList<scalar>> srcFaceWeights_;

                //- For each source face, the coupled target weights
                List<DynamicList<scalar>> tgtFaceWeights_;

                //- Distribution map from source to target faces
                autoPtr<distributionMap> srcFacesMapPtr_;

                //- Distribution map from target to source faces
                autoPtr<distributionMap> tgtFacesMapPtr_;

                //- Map of local source face indices [proc, index]
                autoPtr<List<remote>> localSrcFacesProcPtr_;

                //- Map of local target face indices [proc, index]
                autoPtr<List<remote>> localTgtFacesProcPtr_;


            // Point mapping

                //- Local source to target points
                List<DynamicList<label>> localSrcPointsToTgt_;

                //- Local target to source points
                List<DynamicList<label>> localTgtPointsToSrc_;

                //- For each source point, the coupled target weights
                List<DynamicList<scalar>> srcPointWeights_;

                //- For each source point, the coupled target weights
                List<DynamicList<scalar>> tgtPointWeights_;

                //- Distribution map from source to target points
                autoPtr<distributionMap> srcPointsMapPtr_;

                //- Distribution map from target to source points
                autoPtr<distributionMap> tgtPointsMapPtr_;

                //- Map of local source point indices [proc, index]
                autoPtr<List<remote>> localSrcPointsProcPtr_;

                //- Map of local target point indices [proc, index]
                autoPtr<List<remote>> localTgtPointsProcPtr_;


    // Protected Member Functions

        //- Disallow default bitwise copy construct
        patchToPatchMapping(const patchToPatchMapping&);

        //- Disallow default bitwise assignment
        void operator=(const patchToPatchMapping&);


        // Bound boxes

            //- Create bound box for a single face
            virtual treeBoundBox makeBb
            (
                const face& f,
                const pointField& pts
            ) const;

            //- Create bound box for a single face with normals
            virtual treeBoundBox makeBb
            (
                const face& f,
                const pointField& pts,
                const vectorField& pointNormals
            ) const;

            //- Create bound box for a face list with old points
            treeBoundBox makeBb
            (
                const primitivePatch& patch,
                const pointField& pts0
            ) const;

            //- Create bound box for a face list with old points
            treeBoundBox makeBb
            (
                const primitivePatch& patch,
                const pointField& pts0,
                const label facei
            ) const;

            //- Create bound box for a face list with old points and normals
            treeBoundBox makeBb
            (
                const primitivePatch& patch,
                const pointField& pts0,
                const vectorField& pointNormals,
                const vectorField& pointNormals0
            ) const;

            //- Create bound box for a face list with old points and normals
            treeBoundBox makeBb
            (
                const primitivePatch& patch,
                const pointField& pts0,
                const vectorField& pointNormals,
                const vectorField& pointNormals0,
                const label facei
            ) const;


        // Face intersections

            //- Intersect 2 faces with old points
            virtual bool intersectFaces
            (
                const primitivePatch& srcPatch,
                const pointField& srcPts0,
                const primitivePatch& tgtPatch,
                const pointField& tgtPts0,
                const vectorField& pn,
                const vectorField& pn0,
                const label faceAi,
                const label faceBi
            ) = 0;

            //- Intersect two faces
            bool findOrIntersectFaces
            (
                const primitivePatch& srcPatch,
                const pointField& srcPts0,
                const primitivePatch& tgtPatch,
                const pointField& tgtPts0,
                const vectorField& pointNormals,
                const vectorField& pointNormals0,
                const label faceAi,
                const label faceBi
            );

            //- Intersect a queue of source-target face pairs. Update
            //  completion lists and form a new queue of target-source face
            //  pairs.
            label intersectPatchQueue
            (
                const primitivePatch& srcPatch,
                const pointField& srcPts0,
                const primitivePatch& tgtPatch,
                const pointField& tgtPts0,
                const vectorField& pointNormals,
                const vectorField& pointNormals0,
                const bool isA,
                const DynamicList<labelPair>& queue,
                labelList& faceComplete,
                DynamicList<labelPair>& otherQueue,
                const labelList& otherFaceComplete,
                boolList& otherFaceQueued,
                boolList& otherFaceVisited
            );

            //- Intersect the patches
            void intersectPatches
            (
                const primitivePatch& srcPatch,
                const pointField& pts0A,
                const primitivePatch& tgtPatch,
                const pointField& pts0B,
                const vectorField& pointNormals,
                const vectorField& pointNormals0
            );

        // Parallel operations

            //- Collect faces to send
            labelListList sendTgtPatch
            (
                const primitivePatch& srcPatch,
                const pointField& pts0A,
                const primitivePatch& tgtPatch,
                const pointField& pts0B,
                const vectorField& pointNormals,
                const vectorField& pointNormals0
            ) const;

            //- Collect points and faces to send
            void sendTgtPatch
            (
                const primitivePatch& srcPatch,
                const pointField& pts0A,
                const primitivePatch& tgtPatch,
                const pointField& pts0B,
                const vectorField& pointNormals,
                const vectorField& pointNormals0,
                labelListList& sendPoints,
                labelListList& sendFaces
            ) const;

            //- Create local patch by combining all valid processors
            List<remote> distributePatch
            (
                const distributionMap& map,
                const primitivePatch& patch,
                const pointField& pts0,
                autoPtr<standAlonePatch>& localPatchPtr,
                autoPtr<pointField>& localPoints0Ptr
            );

            //- Create local patch by combining all valid processors
            void distributePatch
            (
                const distributionMap& map,
                const primitivePatch& patch,
                const pointField& pts0,
                List<remote>& localProcPoints,
                List<remote>& localProcFaces,
                autoPtr<standAlonePatch>& localPatchPtr,
                autoPtr<pointField>& localPoints0Ptr
            );


        // Hooks

            //- Initialisation
            virtual void initialise
            (
                const primitivePatch& srcPatch,
                const pointField& srcPts0,
                const primitivePatch& tgtPatch,
                const pointField& tgtPts0,
                const vectorField& pointNormals,
                const vectorField& pointNormals0
            );

            //- Finalise the intersection locally. Trims the local target patch
            //  so that only parts that actually intersect the source remain.
            //  Returns new-to-old map for trimming target-associated data.
            virtual void finaliseLocal
            (
                const primitivePatch& srcPatch,
                const pointField& srcPts0,
                const primitivePatch& tgtPatch,
                const pointField& tgtPts0,
                const vectorField& pointNormals,
                const vectorField& pointNormals0
            );
            virtual labelList finaliseLocalPoints
            (
                const primitivePatch& srcPatch,
                const pointField& srcPts0,
                const primitivePatch& tgtPatch,
                const pointField& tgtPts0,
                const vectorField& pointNormals,
                const vectorField& pointNormals0
            );
            virtual labelList finaliseLocalFaces
            (
                const primitivePatch& srcPatch,
                const pointField& srcPts0,
                const primitivePatch& tgtPatch,
                const pointField& tgtPts0,
                const vectorField& pointNormals,
                const vectorField& pointNormals0
            );

            //- Distribute the source patch so that any information needed by
            //  the target is locally available
            virtual void distributeSrc
            (
                const primitivePatch& patch,
                const pointField& pts0
            );

            //- Send data that resulted from an intersection between the source
            //  patch and a distributed source-local-target patch back to the
            //  original target processes
            virtual void rDistributeTgt
            (
                const primitivePatch& patch,
                const pointField& pts0
            );

            //- Finalising
            virtual labelPair finalise
            (
                const primitivePatch& srcPatch,
                const pointField& srcPts0,
                const primitivePatch& tgtPatch,
                const pointField& tgtPts0,
                const vectorField& pointNormals,
                const vectorField& pointNormals0,
                const transformer& tgtToSrc
            );
            virtual label finalisePoints
            (
                const primitivePatch& srcPatch,
                const pointField& srcPts0,
                const primitivePatch& tgtPatch,
                const pointField& tgtPts0,
                const vectorField& pointNormals,
                const vectorField& pointNormals0,
                const transformer& tgtToSrc
            );
            virtual label finaliseFaces
            (
                const primitivePatch& srcPatch,
                const pointField& srcPts0,
                const primitivePatch& tgtPatch,
                const pointField& tgtPts0,
                const vectorField& pointNormals,
                const vectorField& pointNormals0,
                const transformer& tgtToSrc
            );


public:

    //- Runtime type information
    TypeName("patchToPatchMapping");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            patchToPatchMapping,
            dictionary,
            (
                const primitivePatch& srcPatch,
                const primitivePatch& tgtPatch,
                const dictionary& dict,
                const bool needPoints,
                const bool reverse
            ),
            (srcPatch, tgtPatch, dict, needPoints, reverse)
        );


    // Constructors

        //- Construct from components
        patchToPatchMapping
        (
            const primitivePatch& srcPatch,
            const primitivePatch& tgtPatch,
            const dictionary& dict,
            const bool needPoints,
            const bool reverse
        );


    // Selectors

        //- Select constructed from fluid and solid meshes
        static autoPtr<patchToPatchMapping> New
        (
            const primitivePatch& srcPatch,
            const primitivePatch& tgtPatch,
            const dictionary& dict,
            const bool needPoints,
            const bool reverse
        );


    // Destructor

        virtual ~patchToPatchMapping();


public:


    // Member Functions

        // Access

            //- Return const access to the settings dict
            inline const dictionary& dict() const
            {
                return dict_;
            }

            //- Flag to indicate that the two patches are co-directional and
            //  that the orientation of one should therefore be reversed
            inline bool reverse() const
            {
                return reverse_;
            }

            //- Index of the processor holding all faces of the patchToPatch,
            //  or -1 if spread across multiple processors
            inline label singleProcess() const
            {
                return singleProcess_;
            }

            //- Is this intersection on a single process?
            inline bool isSingleProcess() const
            {
                return singleProcess_ >= 0;
            }

            //- Are point weights used
            bool needPoints() const
            {
                return needPoints_;
            }

            //- Set if point weights are needed
            bool& needPoints()
            {
                return needPoints_;
            }

            //- Return unmapped faces for source patch
            labelList unmappedSrcFaces() const;

            //- Return unmapped faces for target
            labelList unmappedTgtFaces() const;

            //- Return unmapped faces for the given patch
            labelList unmappedFaces(const primitivePatch& patch) const;

            //- Return unmapped points for source patch
            labelList unmappedSrcPoints() const;

            //- Return unmapped points for target
            labelList unmappedTgtPoints() const;

            //- Return unmapped points for the given patch
            labelList unmappedPoints(const primitivePatch& patch) const;


        // Update

            //- Update addressing and weights for the given patches
            void update
            (
                const pointField& srcPts0,
                const pointField& tgtPts0,
                const vectorField& pointNormals,
                const vectorField& pointNormals0,
                const transformer& tgtToSrc = NullObjectRef<transformer>()
            );

            //- Update addressing and weights for the given patches
            void update
            (
                const vectorField& pointNormals,
                const transformer& tgtToSrc = NullObjectRef<transformer>()
            );


            //- General functions that call the virtual implementations

                //- Transfer face field from source to target patch
                template<class Type>
                tmp<Field<Type>> transferFacesToTgt
                (
                    const Field<Type>& srcField,
                    const Field<Type>& unmapped =
                        NullObjectRef<Field<Type>>()
                ) const;
                template<class Type>
                tmp<Field<Type>> transferFacesToTgt
                (
                    const tmp<Field<Type>>& tsrcField,
                    const Field<Type>& unmapped =
                        NullObjectRef<Field<Type>>()
                ) const;

                //- Transfer face field from target to source patch
                template<class Type>
                tmp<Field<Type>> transferFacesToSrc
                (
                    const Field<Type>& tgtField,
                    const Field<Type>& unmapped =
                        NullObjectRef<Field<Type>>()
                ) const;
                template<class Type>
                tmp<Field<Type>> transferFacesToSrc
                (
                    const tmp<Field<Type>>& ttgtField,
                    const Field<Type>& unmapped =
                        NullObjectRef<Field<Type>>()
                ) const;

                //- Transfer point field from source to target patch
                template<class Type>
                tmp<Field<Type>> transferPointsToTgt
                (
                    const Field<Type>& srcField,
                    const Field<Type>& unmapped =
                        NullObjectRef<Field<Type>>()
                ) const;
                template<class Type>
                tmp<Field<Type>> transferPointsToTgt
                (
                    const tmp<Field<Type>>& tsrcField,
                    const Field<Type>& unmapped =
                        NullObjectRef<Field<Type>>()
                ) const;

                //- Transfer point field from target to source patch
                template<class Type>
                tmp<Field<Type>> transferPointsToSrc
                (
                    const Field<Type>& tgtField,
                    const Field<Type>& unmapped =
                        NullObjectRef<Field<Type>>()
                ) const;
                template<class Type>
                tmp<Field<Type>> transferPointsToSrc
                (
                    const tmp<Field<Type>>& ttgtField,
                    const Field<Type>& unmapped =
                        NullObjectRef<Field<Type>>()
                ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "patchToPatchMappingTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
