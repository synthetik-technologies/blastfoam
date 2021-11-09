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
    Foam::immersedMovingWallVelocityFvPatchVectorField

Description


SourceFiles
    immersedMovingWallVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef immersedTemperatureBaseFvPatchScalarField_H
#define immersedTemperatureBaseFvPatchScalarField_H

#include "immersedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class immersedTemperatureBaseFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class immersedTemperatureBaseFvPatchScalarField
:
    public immersedValueFvPatchField<scalar>
{
public:

    //- Runtime type information
    TypeName("immersedMovingWallVelocity");

    // Constructors

        //- Construct from patch and internal field
        immersedTemperatureBaseFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        immersedTemperatureBaseFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given
        // immersedTemperatureBaseFvPatchScalarField onto a new patch
        immersedTemperatureBaseFvPatchScalarField
        (
            const immersedTemperatureBaseFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        immersedTemperatureBaseFvPatchScalarField
        (
            const immersedTemperatureBaseFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        immersedTemperatureBaseFvPatchScalarField
        (
            const immersedTemperatureBaseFvPatchScalarField&,
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
                new immersedTemperatureBaseFvPatchScalarField(*this, iF)
            );
        }


    // Evaluation Functions

        //- Update the boundary values and set the patchInternalField
        virtual void updateCoeffs();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
