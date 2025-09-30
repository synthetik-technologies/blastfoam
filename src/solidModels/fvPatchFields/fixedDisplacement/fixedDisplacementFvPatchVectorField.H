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
    Foam::fixedDisplacementFvPatchVectorField

Description
    Fixed displacement boundary condition with boundary non-orthogonal
    correction for solid solvers.
    Essentially fixedValue with non-orthogonal correction.
    Optionally time-varying.

SourceFiles
    fixedDisplacementFvPatchVectorField.C

Author
    Philip Cardiff, UCD/UT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef fixedDisplacementFvPatchVectorField_H
#define fixedDisplacementFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "Function1.H"
#include "PrimitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class fixedDisplacementFvPatch Declaration
\*---------------------------------------------------------------------------*/

class fixedDisplacementFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
protected:
    // Private Data

        //- Constant total displacement field
        vectorField totalDisp_;

        //- Displacement time series
        autoPtr<Function1<vector>> dispSeries_;

        //- Displacment name
        word displacementName_;

        //- Displacment region
        word regionName_;

        //- Patch interpolator to interpolate face values to the points
        typedef PrimitivePatchInterpolation<primitivePatch>
            primitivePatchInterpolation;
        mutable autoPtr<primitivePatchInterpolation> interpPtr_;

    // Private Functions

        //- Make the patch interpolator
        void makeInterp() const;

        //- Return the patch interpolator
        const primitivePatchInterpolation& interp() const;

        //- Set the point displacement field boundary condition
        void setPointDisplacement(const vectorField& faceDisp);

public:

    //- Runtime type information
    TypeName("fixedDisplacement");


    // Constructors

        //- Construct from patch and internal field
        fixedDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedDisplacementFvPatchVectorField
        //  onto a new patch
        fixedDisplacementFvPatchVectorField
        (
            const fixedDisplacementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fieldMapper&
        );

        //- Construct as copy
        fixedDisplacementFvPatchVectorField
        (
            const fixedDisplacementFvPatchVectorField&
        ) = delete;

        //- Construct as copy setting internal field reference
        fixedDisplacementFvPatchVectorField
        (
            const fixedDisplacementFvPatchVectorField&,
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
                new fixedDisplacementFvPatchVectorField(*this, iF)
            );
        }


    // Destructors

        //- Construct from patch and internal field
        virtual ~fixedDisplacementFvPatchVectorField();


    // Member functions

        // Access functions

        //- Update the coefficients associated with the patch field
            virtual vectorField& totalDisp()
            {
                return totalDisp_;
            }

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Return patch-normal gradient
            //- with non-orthogonal correction regardless of whether snGrad
            //- is corrected or not
            virtual tmp<Field<vector> > snGrad() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<vector> > gradientBoundaryCoeffs() const;

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

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
