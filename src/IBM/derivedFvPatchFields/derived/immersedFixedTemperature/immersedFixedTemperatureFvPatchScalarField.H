/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technology
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    Foam::immersedFixedTemperatureFvPatchScalarField

Description


SourceFiles
    immersedFixedTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef immersedFixedTemperatureFvPatchScalarField_H
#define immersedFixedTemperatureFvPatchScalarField_H

#include "immersedTemperatureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class immersedFixedTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class immersedFixedTemperatureFvPatchScalarField
:
    public immersedTemperatureFvPatchScalarField
{
public:

    //- Runtime type information
    TypeName("immersedFixedTemperature");

    // Constructors

        //- Construct from patch and internal field
        immersedFixedTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        immersedFixedTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given
        // immersedFixedTemperatureFvPatchScalarField onto a new patch
        immersedFixedTemperatureFvPatchScalarField
        (
            const immersedFixedTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        immersedFixedTemperatureFvPatchScalarField
        (
            const immersedFixedTemperatureFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        immersedFixedTemperatureFvPatchScalarField
        (
            const immersedFixedTemperatureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new immersedFixedTemperatureFvPatchScalarField(*this, iF)
            );
        }

    // Member Functions

        // Evaluation functions

            //- Update the boundary values and set the patchInternal field
            virtual void updateCoeffs();

            //- Set the internal values
            virtual void setInternal();


    //- Write
    virtual void write(Ostream&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
