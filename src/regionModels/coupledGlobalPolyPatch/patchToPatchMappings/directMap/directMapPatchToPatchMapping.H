/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    directMapPatchToPatchMapping

Description
    patchToPatchMapping wrapper for direct mapping i.e. conformal
    meshes

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    directMapPatchToPatchMapping.C

\*---------------------------------------------------------------------------*/

#ifndef directMapPatchToPatchMapping_H
#define directMapPatchToPatchMapping_H

#include "patchToPatchMapping.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace patchToPatchMappings
{

/*---------------------------------------------------------------------------*\
                        Class directMapPatchToPatchMapping Declaration
\*---------------------------------------------------------------------------*/

class directMapPatchToPatchMapping
:
    public patchToPatchMapping
{
    // Private data

        //- Face map from zoneB to zoneA
        //  For each zoneB face, give the index of the corresponding zoneA face
        mutable autoPtr<labelIOList> zoneAToZoneBFaceMapPtr_;

        //- Face map from zoneA to zoneB
        //  For each zoneA face, give the index of the corresponding zoneB face
        mutable autoPtr<labelIOList> zoneBToZoneAFaceMapPtr_;

        //- Point map from zoneB to zoneA
        //  For each zoneB point, give the index of the corresponding zoneA
        //  point
        mutable autoPtr<labelIOList> zoneAToZoneBPointMapPtr_;

        //- Point map from zoneA to zoneB
        //  For each zoneA point, give the index of the corresponding zoneB
        //  point
        mutable autoPtr<labelIOList> zoneBToZoneAPointMapPtr_;

        //- List of min edge lengths for interfaces points at fluid side
        mutable autoPtr<scalarField> minEdgeLengthsPtr_;

    // Private static data

        //  Relative tolerance (fraction of minimum edge length used when
        //  checking for matched faces/points
        static const scalar relTol_;

    // Private Member Functions

        //- Check zone face and points sizes are conformal
        void checkZoneSizes() const;

        //- Calculate zoneA to zoneB face map
        void calcZoneAToZoneBFaceMap() const;

        //- Return a reference to the zoneA to zoneB face map
        const labelIOList& zoneAToZoneBFaceMap() const;

        //- Calculate zoneB to zoneA face map
        void calcZoneBToZoneAFaceMap() const;

        //- Return a reference to the zoneB to zoneA face map
        const labelIOList& zoneBToZoneAFaceMap() const;

        //- Calculate zoneA to zoneB point map
        void calcZoneAToZoneBPointMap() const;

        //- Return a reference to the zoneA to zoneB point map
        const labelIOList& zoneAToZoneBPointMap() const;

        //- Calculate zoneB to zoneA point map
        void calcZoneBToZoneAPointMap() const;

        //- Return a reference to the zoneB to zoneA point map
        const labelIOList& zoneBToZoneAPointMap() const;

        //- Calculate minimal edge lengths
        void calcMinEdgeLengths() const;

        //- Return minimal edge lengths
        const scalarField& minEdgeLengths() const;

        //- Transfer/map/interpolate from one zone faces to another zone
        //  faces for Type
        template<class Type>
        void transferFaces
        (
            const standAlonePatch& fromZone,    // from zone
            const standAlonePatch& toZone,      // to zone
            const Field<Type>& fromField,       // from field
            Field<Type>& toField               // to field
        ) const;

        //- Transfer/map/interpolate from one zone points to another zone
        //  points for Type
        template<class Type>
        void transferPoints
        (
            const standAlonePatch& fromZone,    // from zone
            const standAlonePatch& toZone,      // to zone
            const Field<Type>& fromField,       // from field
            Field<Type>& toField               // to field
        ) const;

        //- Disallow default bitwise copy construct
        directMapPatchToPatchMapping
        (
            const directMapPatchToPatchMapping&
        );

        //- Disallow default bitwise assignment
        void operator=(const directMapPatchToPatchMapping&);


public:

    //- Runtime type information
    TypeName("directMap");


    // Constructors

        //- Construct from components
        directMapPatchToPatchMapping
        (
            const dictionary& dict,
            const primitivePatch& patchA,
            const primitivePatch& patchB,
            const globalPolyPatch& globalPatchA,
            const globalPolyPatch& globalPatchB
        );


    // Destructor

        virtual ~directMapPatchToPatchMapping()
        {}


    // Member Functions

        // Edit

            //- Transfer/map/interpolate from one zone faces to another zone
            //  faces for scalars
            virtual void transferFaces
            (
                const standAlonePatch& fromPatch,      // from zone
                const standAlonePatch& toPatch,        // to zone
                const Field<scalar>& fromField,       // from field
                Field<scalar>& toField                // to field
            ) const;

            //- Transfer/map/interpolate from one zone faces to another zone
            //  faces for vectors
            virtual void transferFaces
            (
                const standAlonePatch& fromPatch,      // from zone
                const standAlonePatch& toPatch,        // to zone
                const Field<vector>& fromField,       // from field
                Field<vector>& toField                // to field
            ) const;

            //- Transfer/map/interpolate from one zone faces to another zone
            //  faces for symmTensor
            virtual void transferFaces
            (
                const standAlonePatch& fromPatch,      // from zone
                const standAlonePatch& toPatch,        // to zone
                const Field<symmTensor>& fromField,   // from field
                Field<symmTensor>& toField            // to field
            ) const;

            //- Transfer/map/interpolate from one zone faces to another zone
            //  faces for sphericalTensor
            virtual void transferFaces
            (
                const standAlonePatch& fromPatch,      // from zone
                const standAlonePatch& toPatch,        // to zone
                const Field<sphericalTensor>& fromField, // from field
                Field<sphericalTensor>& toField       // to field
            ) const;

            //- Transfer/map/interpolate from one zone faces to another zone
            //  faces for tensor
            virtual void transferFaces
            (
                const standAlonePatch& fromPatch,      // from zone
                const standAlonePatch& toPatch,        // to zone
                const Field<tensor>& fromField,       // from field
                Field<tensor>& toField                // to field
            ) const;


            //- Transfer/map/interpolate from one zone points to another zone
            //  points for scalars
            virtual void transferPoints
            (
                const standAlonePatch& fromPatch,      // from zone
                const standAlonePatch& toPatch,        // to zone
                const Field<scalar>& fromField,       // from field
                Field<scalar>& toField                // to field
            ) const;

            //- Transfer/map/interpolate from one zone points to another zone
            //  points for vectors
            virtual void transferPoints
            (
                const standAlonePatch& fromPatch,      // from zone
                const standAlonePatch& toPatch,        // to zone
                const Field<vector>& fromField,       // from field
                Field<vector>& toField                // to field
            ) const;

            //- Transfer/map/interpolate from one zone points to another zone
            //  points for symmTensors
            virtual void transferPoints
            (
                const standAlonePatch& fromPatch,      // from zone
                const standAlonePatch& toPatch,        // to zone
                const Field<symmTensor>& fromField,   // from field
                Field<symmTensor>& toField            // to field
            ) const;

            //- Transfer/map/interpolate from one zone points to another zone
            //  points for sphericalTensor
            virtual void transferPoints
            (
                const standAlonePatch& fromPatch,      // from zone
                const standAlonePatch& toPatch,        // to zone
                const Field<sphericalTensor>& fromField, // from field
                Field<sphericalTensor>& toField       // to field
            ) const;

            //- Transfer/map/interpolate from one zone points to another zone
            //  points for tensors
            virtual void transferPoints
            (
                const standAlonePatch& fromPatch,      // from zone
                const standAlonePatch& toPatch,        // to zone
                const Field<tensor>& fromField,       // from field
                Field<tensor>& toField                // to field
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchToPatchMappings
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "directMapPatchToPatchMappingTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
