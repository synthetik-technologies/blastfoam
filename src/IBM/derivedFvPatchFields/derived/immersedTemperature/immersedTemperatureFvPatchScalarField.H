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
    Foam::immersedTemperatureFvPatchScalarField

Description


SourceFiles
    immersedTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef immersedTemperatureFvPatchScalarField_H
#define immersedTemperatureFvPatchScalarField_H

#include "immersedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class basicThermo;

/*---------------------------------------------------------------------------*\
            Class immersedTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class immersedTemperatureFvPatchScalarField
:
    public immersedValueFvPatchField<scalar>
{
protected:
// Protected member data

    //- Is the internal field energy
    bool isEnergy_;


    //- Set the patch energy
    void setPatchHE();

    //- Set the internal energy
    void setInternalHE();


public:

    //- Runtime type information
//     TypeName("immersedTemperature");

    // Constructors

        //- Construct from patch and internal field
        immersedTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        immersedTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );

        //- Construct by mapping the given
        // immersedTemperatureFvPatchScalarField onto a new patch
        immersedTemperatureFvPatchScalarField
        (
            const immersedTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        immersedTemperatureFvPatchScalarField
        (
            const immersedTemperatureFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        immersedTemperatureFvPatchScalarField
        (
            const immersedTemperatureFvPatchScalarField&,
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
                new immersedTemperatureFvPatchScalarField(*this, iF)
            );
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
