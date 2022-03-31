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
    analyticalPlateHoleTractionFvPatchVectorField

Description
    Traction boundary condition for plateHole test case to enforce
    analytical traction at the boundaries

Author
   Philip Cardiff, UCD.
   Zeljko Tukovic, Zagreb.

SourceFiles
    analyticalPlateHoleTractionFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef analyticalPlateHoleTractionFvPatchVectorField_H
#define analyticalPlateHoleTractionFvPatchVectorField_H

#include "fvPatchFields.H"
#include "solidTractionFvPatchVectorField.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class analyticalPlateHoleTractionFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class analyticalPlateHoleTractionFvPatchVectorField
:
    public solidTractionFvPatchVectorField
{
private:

    // Private Data

        //- Far field traction
        //  It is assumed that it acts in the positive X direction
        const scalar T_;

        //- Hole radius
        const scalar holeR_;


    // Private Member Functions

        //- Analytical solution
        symmTensor plateHoleSolution(const vector& C);


public:

    //- Runtime type information
    TypeName("analyticalPlateHoleTraction");


    // Constructors

        //- Construct from patch and internal field
        analyticalPlateHoleTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        analyticalPlateHoleTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  analyticalPlateHoleTractionFvPatchVectorField onto a new patch
        analyticalPlateHoleTractionFvPatchVectorField
        (
            const analyticalPlateHoleTractionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        analyticalPlateHoleTractionFvPatchVectorField
        (
            const analyticalPlateHoleTractionFvPatchVectorField&
        ) = delete;

        //- Construct as copy setting internal field reference
        analyticalPlateHoleTractionFvPatchVectorField
        (
            const analyticalPlateHoleTractionFvPatchVectorField&,
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
                new analyticalPlateHoleTractionFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

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

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
