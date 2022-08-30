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
    rbfPatchToPatchMapping

Description
    patchToPatchMapping wrapper using radial basis functions

Author
    Philip Cardiff, UCD. All rights reserved.
    This class is a wrapper for the code from David Blom

SourceFiles
    rbfPatchToPatchMapping.C

\*---------------------------------------------------------------------------*/

#ifndef rbfPatchToPatchMapping_H
#define rbfPatchToPatchMapping_H

#include "patchToPatchMapping.H"
#include "RBFInterpolation.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace patchToPatchMappings
{

/*---------------------------------------------------------------------------*\
                        Class rbfPatchToPatchMapping Declaration
\*---------------------------------------------------------------------------*/

class rbfPatchToPatchMapping
:
    public patchToPatchMapping
{
    // Private data

        //- Dictionary
        const dictionary& dict_;

        //- zoneA to zoneB interpolator
        mutable autoPtr<RBFInterpolation> zoneAToZoneBInterpolatorPtr_;

        //- List of solid zone to fluid zone interpolators
        mutable autoPtr<RBFInterpolation> zoneBToZoneAInterpolatorPtr_;

    // Private Member Functions

        //- Make zoneA to zoneB interpolator
        void makeZoneAToZoneBInterpolator() const;

        //- Return reference to zoneA to zoneB interpolator
        const RBFInterpolation& zoneAToZoneBInterpolator() const;

        //- Make zoneB to zoneA interpolator
        void makeZoneBToZoneAInterpolator() const;

        //- Return reference to zoneB to zoneA interpolator
        const RBFInterpolation& zoneBToZoneAInterpolator() const;

        //- Transfer/map/interpolate from one zone faces to another zone
        //  faces for Type
        template<class Type>
        void transferFaces
        (
            const standAlonePatch& fromZone,    // from zone
            const standAlonePatch& toZone,      // to zone
            const Field<Type>& fromField,       // from field
            Field<Type>& toField                // to field
        ) const;

        //- Transfer/map/interpolate from one zone points to another zone
        //  points for Type
        template<class Type>
        void transferPoints
        (
            const standAlonePatch& fromZone,    // from zone
            const standAlonePatch& toZone,      // to zone
            const Field<Type>& fromField,       // from field
            Field<Type>& toField                // to field
        ) const;

        //- Disallow default bitwise copy construct
        rbfPatchToPatchMapping
        (
            const rbfPatchToPatchMapping&
        );

        //- Disallow default bitwise assignment
        void operator=(const rbfPatchToPatchMapping&);


public:

    //- Runtime type information
    TypeName("RBF");


    // Constructors

        //- Construct from components
        rbfPatchToPatchMapping
        (
            const dictionary& dict,
            const primitivePatch& patchA,
            const primitivePatch& patchB,
            const globalPolyPatch& globalPatchA,
            const globalPolyPatch& globalPatchB
        );


    // Destructor

        virtual ~rbfPatchToPatchMapping()
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

#ifdef NoRepository
#   include "rbfPatchToPatchMappingTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
