/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    cohesiveZoneModel

Description
    Cohesive zone model class stores the cohesive properties, such as strengths
    and toughness, and has two main functions:

    (1) cohesiveTractions
        Return the tractions to apply to the crack patch, which may be a
        function of the crack patch displacements, energies, etc.

    (2) updateEnergies
        Calculate the current dissipated fracture energies;
Authors
    Philip Cardif, UCD/UT. All rights reserved.

SourceFiles
    cohesiveZoneModel.C
    newCohesiveZoneModel.C

\*---------------------------------------------------------------------------*/

#ifndef cohesiveZoneModel_H
#define cohesiveZoneModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "tmp.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Class forward declaration
    class solidCohesiveFvPatchVectorField;

/*---------------------------------------------------------------------------*\
                         Class cohesiveZoneModel Declaration
\*---------------------------------------------------------------------------*/

class cohesiveZoneModel
{
    // Private data

        //- Name
        const word name_;

        //- Const reference to the fvPatch
        const fvPatch& patch_;

        //- Traction field on the mesh
        mutable surfaceVectorField traction_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const cohesiveZoneModel&);


public:

    //- Runtime type information
    TypeName("cohesiveZoneModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            cohesiveZoneModel,
            dictionary,
            (
                const word name,
                const fvPatch& patch,
                const dictionary& dict
            ),
            (name, patch, dict)
        );


    // Selectors

        //- Return a reference to the selected meshFailure model
        static autoPtr<cohesiveZoneModel> New
        (
            const word& name,
            const fvPatch& patch,
            const dictionary& dict
        );


    // Constructors

        //- Construct from dictionary
        cohesiveZoneModel
        (
            const word& name,
            const fvPatch& patch,
            const dictionary& dict
        );

        //- Construct as copy
        cohesiveZoneModel(const cohesiveZoneModel&);

        //- Construct and return a clone
        virtual autoPtr<cohesiveZoneModel> clone() const = 0;


    // Destructor

        virtual ~cohesiveZoneModel()
        {}


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper& mapper
            ) = 0;

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const solidCohesiveFvPatchVectorField& sc,
                const labelList& addr
            ) = 0;

        // Access functions

            //- Return name
            const word& name() const
            {
                return name_;
            }

            //- Return const reference to the fvPatch
            const fvPatch& patch() const
            {
                return patch_;
            }

            //- Return const reference to the fvMesh
            const fvMesh& mesh() const
            {
                return patch_.boundaryMesh().mesh();
            }

            //- Return a const reference to the current mesh traction field
            const surfaceVectorField& meshTraction() const;

            //- Update the mesh traction field
            void updateMeshTraction() const;

            //- Update old values
            virtual void updateOldFields() = 0;

            //- Update fracture energy
            virtual void updateEnergy
            (
                const vectorField& traction,
                const vectorField& delta
            ) = 0;

            //- Update cohesive traction field
            virtual void updateTraction
            (
                vectorField& traction,
                const vectorField& delta
            ) = 0;

            //- Return traction fraction for undamaged faces
            //  This function is called by the faceBreakerLaw to determine which
            //  internal faces should be add to the cohesive zone
            virtual tmp<surfaceScalarField>
                initiationTractionFraction() const = 0;

            //- Crack and damage field for visualisation
            //  Cracked faces set to 2, damaged faces set to 1
            virtual tmp<scalarField> crackingAndDamage() const = 0;

            //- Write dictionary
            virtual void write(Ostream&) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
