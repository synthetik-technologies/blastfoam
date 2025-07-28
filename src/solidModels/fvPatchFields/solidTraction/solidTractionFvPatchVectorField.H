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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    solidTractionFvPatchVectorField

Description
    Fixed traction boundary condition for solid solvers.

Author
    Philip Cardiff, UCD.
    Zeljko Tukovic, FSB Zagreb.

SourceFiles
    solidTractionFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef solidTractionFvPatchVectorField_H
#define solidTractionFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "tractionBase.H"
#include "Switch.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class solidTractionFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class solidTractionFvPatchVectorField
:
    public fixedGradientFvPatchVectorField,
    public tractionBase
{
private:
    // Private Data

        //- Traction time series
        autoPtr<Function1<vector>> tractionSeries_;

        //- Pressure time series
        autoPtr<Function1<scalar>> pressureSeries_;

        //- Second order accuracy correction
        Switch secondOrder_;

        //- Optional: under-relaxation factor
        autoPtr<Function1<scalar>> relaxFac_;


public:

    //- Runtime type information
    TypeName("solidTraction");


    // Constructors

        //- Construct from patch and internal field
        solidTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&,
            const bool mustRead = true
        );

        //- Construct by mapping given
        //  solidTractionFvPatchVectorField onto a new patch
        solidTractionFvPatchVectorField
        (
            const solidTractionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fieldMapper&
        );

        //- Construct as copy
        solidTractionFvPatchVectorField
        (
            const solidTractionFvPatchVectorField&
        ) = delete;

        //- Construct as copy setting internal field reference
        solidTractionFvPatchVectorField
        (
            const solidTractionFvPatchVectorField&,
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
                new solidTractionFvPatchVectorField(*this, iF)
            );
        }



    // Member functions

        // Access

            //- Return if second order
            inline bool secondOrder() const
            {
                return secondOrder_;
            }


        // Mapping functions

            //- Map the given fvPatchField onto this fvPatchField
            virtual void map
            (
                const fvPatchField<vector>&,
                const fieldMapper&
            );

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchField<vector>&);


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Update pressure and traction
        virtual bool updateFields();

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType = Pstream::commsTypes::blocking
        );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
