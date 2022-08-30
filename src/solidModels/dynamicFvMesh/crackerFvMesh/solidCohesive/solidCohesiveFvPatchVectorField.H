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
    solidCohesiveFvPatchVectorField

Description
    Coheisve zone model boundary for arbitrary crack procedure, where the actual
    cohesiveZoneModel is run-time selectable.

Author
    Philip Cardiff, UCD/UT. All rights reserved.
    Zeljko Tukovic, FSB Zaregb. All rights reserved.

SourceFiles
    solidCohesiveFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef solidCohesiveFvPatchVectorField_H
#define solidCohesiveFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "cohesiveFvPatch.H"
#include "Switch.H"
#include "cohesiveZoneModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class solidCohesiveFvPatch Declaration
\*---------------------------------------------------------------------------*/

class solidCohesiveFvPatchVectorField
:
    public fixedGradientFvPatchVectorField
{

    // Private Data

        //- Name of the displacement field
        const word fieldName_;

        //- Run-time selectable cohesive zone model
        mutable autoPtr<cohesiveZoneModel> cohesiveZoneModelPtr_;

        //- Current traction field
        vectorField traction_;

        //- Current pressure field - this is used to represent pore-pressure
        scalarField pressure_;

        //- Current face deltas
        //  The difference in displacement between a face and its shadow face
        vectorField delta_;

        //- Under-relaxation factor for deltas
        const scalar relaxationFactor_;

        //- Under-relaxation factor for tractions
        const scalar relaxationFactorTrac_;

        //- Current time index
        label curTimeIndex_;

        //- Is the solevr using an incremental approach
        //  Is it solving for DD instead of D
        const bool incremental_;

        // Copy of patch dictionary
        const dictionary dict_;


    // Private Member Functions

        //- Update delta
        void updateDelta();

        //- Set the cohesive zone model
        void calcCohesiveZone() const;


public:

    //- Runtime type information
    TypeName("solidCohesive");


    // Constructors

        //- Construct from patch and internal field
        solidCohesiveFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidCohesiveFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  solidCohesiveFvPatchVectorField onto a new patch
        solidCohesiveFvPatchVectorField
        (
            const solidCohesiveFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solidCohesiveFvPatchVectorField
        (
            const solidCohesiveFvPatchVectorField&
        ) = delete;

        //- Construct as copy setting internal field reference
        solidCohesiveFvPatchVectorField
        (
            const solidCohesiveFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new solidCohesiveFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return true if this patch field fixes a value.
            virtual bool fixesValue() const
            {
                return false;
            }

            //- Return const reference to cohesive zone model
            const cohesiveZoneModel& cohesiveZone() const;

            //- Return reference to cohesive zone model
            cohesiveZoneModel& cohesiveZone();

            //- Return current traction
            vectorField& traction()
            {
                return traction_;
            }

            //- Return const reference to pressure
            const scalarField& pressure() const
            {
                return pressure_;
            }

            //- Return reference to pressure
            scalarField& pressure()
            {
                return pressure_;
            }

            //- Return const reference to current deltas
            const vectorField& delta() const
            {
                return delta_;
            }

            //- Return reference to current deltas
            vectorField& delta()
            {
                return delta_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType =
                Pstream::commsTypes::blocking
        );

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchVectorField&) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
