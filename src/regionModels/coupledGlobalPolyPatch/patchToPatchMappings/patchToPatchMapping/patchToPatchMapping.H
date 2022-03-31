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
    patchToPatchMapping

Description
    Virtual base class for mapping/interpolating between two interfaces

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    patchToPatchMapping.C

\*---------------------------------------------------------------------------*/

#ifndef patchToPatchMapping_H
#define patchToPatchMapping_H

#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "primitivePatch.H"
#include "globalPolyPatch.H"

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

        //- Is it required that the patches match
        bool requireMatch_;


private:
    // Private data

        //- Const reference to the patchA
        const primitivePatch& patchA_;

        //- Const reference to the patchB
        const primitivePatch& patchB_;

        //- Const reference to the globalPatchA
        const globalPolyPatch& globalPatchA_;

        //- Const reference to the globalPatchB
        const globalPolyPatch& globalPatchB_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        patchToPatchMapping(const patchToPatchMapping&);

        //- Disallow default bitwise assignment
        void operator=(const patchToPatchMapping&);


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
                const dictionary& dict,
                const primitivePatch& patchA,
                const primitivePatch& patchB,
                const globalPolyPatch& globalPatchA,
                const globalPolyPatch& globalPatchB
            ),
            (dict, patchA, patchB, globalPatchA, globalPatchB)
        );


    // Constructors

        //- Construct from components
        patchToPatchMapping
        (
            const word& type,
            const dictionary& dict,
            const primitivePatch& patchA,
            const primitivePatch& patchB,
            const globalPolyPatch& globalPatchA,
            const globalPolyPatch& globalPatchB
        );


    // Selectors

        //- Select constructed from fluid and solid meshes
        static autoPtr<patchToPatchMapping> New
        (
            const dictionary& dict,
            const primitivePatch& patchA,
            const primitivePatch& patchB,
            const globalPolyPatch& globalPatchA,
            const globalPolyPatch& globalPatchB
        );


    // Destructor

        virtual ~patchToPatchMapping();


public:


    // Member Functions

        // Access

            //- Return const access to the settings dict
            const dictionary& dict() const
            {
                return dict_;
            }

            //- Const access to the patchA
            const primitivePatch& patchA() const
            {
                return patchA_;
            }

            //- Const access to the patchB
            const primitivePatch& patchB() const
            {
                return patchB_;
            };

            //- Const access to the globalPatchA
            const globalPolyPatch& globalPatchA() const
            {
                return globalPatchA_;
            }

            //- Const access to the globalPatchB
            const globalPolyPatch& globalPatchB() const
            {
                return globalPatchB_;
            };

            //- Const access to the primitivePatch of globalPolyPatchA
            const standAlonePatch& zoneA() const
            {
                return globalPatchA_.globalPatch();
            }

            //- Const access to the primitivePatch of globalPolyPatchB
            const standAlonePatch& zoneB() const
            {
                return globalPatchB_.globalPatch();
            };

            //- Check field sizes are correct before interpolation/mapping
            void checkFieldSizes
            (
                const label fromPatchSize,
                const label toPatchSize,
                const label fromFieldSize,
                const label toFieldSize
            ) const;

            //- General functions that call the virtual implementations

                //- Transfer faces from the given patch to the other patch
                template<class Type>
                tmp<Field<Type>> transferFaces
                (
                    const standAlonePatch& fromPatch,
                    const Field<Type>& fromField
                ) const;

                template<class Type>
                tmp<Field<Type>> transferFaces
                (
                    const standAlonePatch& fromPatch,
                    const tmp<Field<Type>>& fromField
                ) const;

                //- Transfer points from the given patch to the other patch
                template<class Type>
                tmp<Field<Type>> transferPoints
                (
                    const standAlonePatch& fromPatch,
                    const Field<Type>& fromField
                ) const;

                template<class Type>
                tmp<Field<Type>> transferPoints
                (
                    const standAlonePatch& fromPatch,
                    const tmp<Field<Type>>& fromField
                ) const;


            //  Note: unfortunatley, a template virtual function is not allowed
            //  in C++, so instead we will make an explicit implementation for
            //  scalars and vectors; derived classes can use templates to make
            //  their lives easier if desired

            //- Transfer/map/interpolate from one zone faces to another zone
            //  faces for scalars
            virtual void transferFaces
            (
                const standAlonePatch& fromPatch,      // from zone
                const standAlonePatch& toPatch,        // to zone
                const Field<scalar>& fromField,       // from field
                Field<scalar>& toField                // to field
            ) const = 0;

            //- Transfer/map/interpolate from one zone faces to another zone
            //  faces for vectors
            virtual void transferFaces
            (
                const standAlonePatch& fromPatch,      // from zone
                const standAlonePatch& toPatch,        // to zone
                const Field<vector>& fromField,       // from field
                Field<vector>& toField                // to field
            ) const = 0;

            //- Transfer/map/interpolate from one zone faces to another zone
            //  faces for symmTensor
            virtual void transferFaces
            (
                const standAlonePatch& fromPatch,      // from zone
                const standAlonePatch& toPatch,        // to zone
                const Field<symmTensor>& fromField,   // from field
                Field<symmTensor>& toField            // to field
            ) const = 0;

            //- Transfer/map/interpolate from one zone faces to another zone
            //  faces for sphericalTensor
            virtual void transferFaces
            (
                const standAlonePatch& fromPatch,      // from zone
                const standAlonePatch& toPatch,        // to zone
                const Field<sphericalTensor>& fromField, // from field
                Field<sphericalTensor>& toField       // to field
            ) const = 0;

            //- Transfer/map/interpolate from one zone faces to another zone
            //  faces for tensor
            virtual void transferFaces
            (
                const standAlonePatch& fromPatch,      // from zone
                const standAlonePatch& toPatch,        // to zone
                const Field<tensor>& fromField,       // from field
                Field<tensor>& toField                // to field
            ) const = 0;


            //- Transfer/map/interpolate from one zone points to another zone
            //  points for scalars
            virtual void transferPoints
            (
                const standAlonePatch& fromPatch,      // from zone
                const standAlonePatch& toPatch,        // to zone
                const Field<scalar>& fromField,       // from field
                Field<scalar>& toField                // to field
            ) const = 0;

            //- Transfer/map/interpolate from one zone points to another zone
            //  points for vectors
            virtual void transferPoints
            (
                const standAlonePatch& fromPatch,      // from zone
                const standAlonePatch& toPatch,        // to zone
                const Field<vector>& fromField,       // from field
                Field<vector>& toField                // to field
            ) const = 0;

            //- Transfer/map/interpolate from one zone points to another zone
            //  points for symmTensors
            virtual void transferPoints
            (
                const standAlonePatch& fromPatch,      // from zone
                const standAlonePatch& toPatch,        // to zone
                const Field<symmTensor>& fromField,   // from field
                Field<symmTensor>& toField            // to field
            ) const = 0;

            //- Transfer/map/interpolate from one zone points to another zone
            //  points for sphericalTensor
            virtual void transferPoints
            (
                const standAlonePatch& fromPatch,      // from zone
                const standAlonePatch& toPatch,        // to zone
                const Field<sphericalTensor>& fromField, // from field
                Field<sphericalTensor>& toField       // to field
            ) const = 0;

            //- Transfer/map/interpolate from one zone points to another zone
            //  points for tensors
            virtual void transferPoints
            (
                const standAlonePatch& fromPatch,      // from zone
                const standAlonePatch& toPatch,        // to zone
                const Field<tensor>& fromField,       // from field
                Field<tensor>& toField                // to field
            ) const = 0;
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
