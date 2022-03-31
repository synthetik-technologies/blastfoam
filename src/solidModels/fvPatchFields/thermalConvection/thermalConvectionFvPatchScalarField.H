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
    thermalConvectionFvPatchScalarField

Description
    thermalConvection boundary condition for temperature
    free/forced convection boundary.

    In progress: non-orthogonal correction not implemented yet.

SourceFiles
    thermalConvectionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/


#ifndef thermalConvectionFvPatchScalarFields_H
#define thermalConvectionFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fixedFluxPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class thermalConvectionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Name of the thermal diffusivity field
        word DTName_;

        //- Heat transfer coefficient
        scalarField alpha_;

        //- Temperature at infinity
        dimensionedScalar Tinf_;

public:

    //- Runtime type information
    TypeName("thermalConvection");


    // Constructors

        //- Construct from patch and internal field
        thermalConvectionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        thermalConvectionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedFluxPressureFvPatchScalarField onto
        //  a new patch
        thermalConvectionFvPatchScalarField
        (
            const thermalConvectionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        thermalConvectionFvPatchScalarField
        (
            const thermalConvectionFvPatchScalarField&
        ) = delete;

        //- Construct as copy setting internal field reference
        thermalConvectionFvPatchScalarField
        (
            const thermalConvectionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new thermalConvectionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // //- Update the coefficients associated with the patch field
        // virtual void updateCoeffs();

        //- Return gradient at boundary
        virtual tmp<Field<scalar> > snGrad() const;

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType =
                Pstream::commsTypes::blocking
        );

        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<Field<scalar> > gradientInternalCoeffs() const;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<Field<scalar> > gradientBoundaryCoeffs() const;

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<scalar>&,
                const labelList&
            );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
