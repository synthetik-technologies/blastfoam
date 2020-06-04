/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::fixedFaceFluxPressureFvPatchScalarField

Description
    This boundary condition sets the pressure gradient to the provided value
    such that the flux on the boundary is that specified by the velocity
    boundary condition.

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fixedFluxPressure;
    }
    \endverbatim

See also
    Foam::fixedGradientFvPatchField

SourceFiles
    fixedFaceFluxPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedFaceFluxPressureFvPatchScalarField_H
#define fixedFaceFluxPressureFvPatchScalarField_H

#include "fvPatchFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class fixedFaceFluxPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fixedFaceFluxPressureFvPatchScalarField
:
    public zeroGradientFvPatchScalarField
{
    // Private data

        label refFace_;
        scalar refP_;


public:

    //- Runtime type information
    TypeName("fixedFaceFluxPressure");


    // Constructors

        //- Construct from patch and internal field
        fixedFaceFluxPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedFaceFluxPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedFaceFluxPressureFvPatchScalarField onto
        //  a new patch
        fixedFaceFluxPressureFvPatchScalarField
        (
            const fixedFaceFluxPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedFaceFluxPressureFvPatchScalarField
        (
            const fixedFaceFluxPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedFaceFluxPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedFaceFluxPressureFvPatchScalarField
        (
            const fixedFaceFluxPressureFvPatchScalarField&,
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
                new fixedFaceFluxPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

         virtual bool fixesValue() const
        {
            return true;
        }

        //- Update the patch pressure gradient field from the given snGradp
//         virtual void updateCoeffs(const scalarField& snGradp);

        //- Update the patch pressure gradient field
//         virtual void updateCoeffs();

        void evaluate(const Pstream::commsTypes);

        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
        virtual tmp<scalarField> valueInternalCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
        virtual tmp<scalarField> valueBoundaryCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
        virtual tmp<scalarField> gradientInternalCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
        virtual tmp<scalarField> gradientBoundaryCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
