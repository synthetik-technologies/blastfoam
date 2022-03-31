/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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
    Foam::convectiveFluxFvPatchScalarField

Description
    Foam::convectiveFluxFvPatchScalarField

SourceFiles
    convectiveFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef convectiveFluxFvPatchScalarFields_H
#define convectiveFluxFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fixedFluxPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class convectiveFluxFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        // //- Name of the velocity field
        // word UName_;

        // //- Name of the flux transporting the field
        // word phiName_;

        // //- Name of the density field used to normalise the mass flux
        // //  if neccessary
        // word rhoName_;

        // //- Is the pressure adjoint, i.e. has the opposite sign
        // Switch adjoint_;

        //- Heat transfer coefficient
        scalarField alpha_;

//         //- Conductivity
//         scalarField lambda_;

        //- Temperature at infinity
        dimensionedScalar Tinf_;

public:

    //- Runtime type information
    TypeName("convectiveFlux");


    // Constructors

        //- Construct from patch and internal field
        convectiveFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        convectiveFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedFluxPressureFvPatchScalarField onto
        //  a new patch
        convectiveFluxFvPatchScalarField
        (
            const convectiveFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        convectiveFluxFvPatchScalarField
        (
            const convectiveFluxFvPatchScalarField&
        ) = delete;

        //- Construct as copy setting internal field reference
        convectiveFluxFvPatchScalarField
        (
            const convectiveFluxFvPatchScalarField&,
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
                new convectiveFluxFvPatchScalarField(*this, iF)
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

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
