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
    Foam::immersedZeroGradientFvPatchField

Description


SourceFiles
    immersedZeroGradientFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef immersedZeroGradientFvPatchField_H
#define immersedZeroGradientFvPatchField_H

#include "volFields.H"
#include "immersedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class immersedZeroGradientFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class immersedZeroGradientFvPatchField
:
    public immersedValueFvPatchField<Type>
{
public:

    //- Runtime type information
    TypeName("immersedZeroGradient");

    // Constructors

        //- Construct from patch and internal field
        immersedZeroGradientFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        immersedZeroGradientFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given
        // immersedZeroGradientFvPatchField<Type> onto a new patch
        immersedZeroGradientFvPatchField
        (
            const immersedZeroGradientFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        immersedZeroGradientFvPatchField
        (
            const immersedZeroGradientFvPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        immersedZeroGradientFvPatchField
        (
            const immersedZeroGradientFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new immersedZeroGradientFvPatchField<Type>(*this, iF)
            );
        }


    // Evaluation functions

        //- Update the boundary values and set the patchInternal field
        virtual void updateCoeffs();

        //- Set the internal values
        virtual void setInternal();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "immersedZeroGradientFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
