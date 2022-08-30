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
    Foam::fixedRotationFvPatchVectorField

Description
    Fixed rotation boundary condition with boundary non-orthogonal
    correction for stress solvers.

    Optional: time-varying angle series.

SourceFiles
    fixedRotationFvPatchVectorField.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef fixedRotationFvPatchVectorField_H
#define fixedRotationFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class fixedRotationFvPatch Declaration
\*---------------------------------------------------------------------------*/

class fixedRotationFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

    // Private Data

        //- Rotation in degrees
        scalar rotationAngle_;

        //- Rotation axis
        const vector rotationAxis_;

        //- Origin of rotation
        vector rotationOrigin_;

        //- Initial face centres
        const vectorField origFaceCentres_;

        //- Initial patch points
        const vectorField origPatchPoints_;

        //- Optional angle time series
        autoPtr<Function1<scalar>> angleSeries_;

        //- Optional translation series
        autoPtr<Function1<vector>> dispSeries_;

        //- Optional origin series
        autoPtr<Function1<vector>> originSeries_;


public:

    //- Runtime type information
    TypeName("fixedRotation");


    // Constructors

        //- Construct from patch and internal field
        fixedRotationFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedRotationFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedRotationFvPatchVectorField
        //  onto a new patch
        fixedRotationFvPatchVectorField
        (
            const fixedRotationFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedRotationFvPatchVectorField
        (
            const fixedRotationFvPatchVectorField&
        ) = delete;

        //- Construct as copy setting internal field reference
        fixedRotationFvPatchVectorField
        (
            const fixedRotationFvPatchVectorField&,
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
                new fixedRotationFvPatchVectorField(*this, iF)
            );
        }



    // Member functions

        // Evaluation functions

            //- Return patch-normal gradient
            //- with non-orthogonal correction regardless of whether snGrad
            //- is corrected or not
            virtual tmp<Field<vector> > snGrad() const;

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<vector> > gradientBoundaryCoeffs() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
