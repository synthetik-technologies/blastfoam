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
#include "solidTractionFvPatchVectorField.H"
#include "cohesiveFvPatch.H"
#include "Switch.H"
#include "cohesiveZoneModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class SolidCohesiveFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

template<class SolidTraction>
class SolidCohesiveFvPatchVectorField
:
    public SolidTraction,
    public cohesiveZoneModelMaster
{

    // Private Data

        //- Current face deltas
        //  The difference in displacement between a face and its shadow face
        vectorField delta_;

        //- Under-relaxation factor for deltas
        const scalar relaxationFactorDelta_;

        //- Under-relaxation factor for tractions
        const scalar relaxationFactorTraction_;

        //- Current time index
        label curTimeIndex_;



    // Private Member Functions

        //- Update delta
        void updateDelta();


public:

    TypeName("solidCohesiveTraction");

    // Constructors

        //- Construct from patch and internal field
        SolidCohesiveFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        SolidCohesiveFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  SolidCohesiveFvPatchVectorField onto a new patch
        SolidCohesiveFvPatchVectorField
        (
            const SolidCohesiveFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fieldMapper&
        );

        //- Construct as copy
        SolidCohesiveFvPatchVectorField
        (
            const SolidCohesiveFvPatchVectorField&
        ) = delete;

        //- Construct as copy setting internal field reference
        SolidCohesiveFvPatchVectorField
        (
            const SolidCohesiveFvPatchVectorField&,
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
                new SolidCohesiveFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return true if this patch field fixes a value.
            virtual bool fixesValue() const
            {
                return false;
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

            //- Map the given fvPatchField onto this fvPatchField
            virtual void map
            (
                const fvPatchVectorField& ptf,
                const fieldMapper&
            );

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset
            (
                const fvPatchVectorField&
            );


        //- Update the coefficients associated with the patch field
        virtual bool updateFields();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchVectorField&) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "solidCohesiveFvPatchVectorField.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
