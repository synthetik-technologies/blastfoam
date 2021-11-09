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
    Foam::immersedBasicSymmetryFvPatchField

Description


SourceFiles
    immersedBasicSymmetryFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBasicSymmetryFvPatchField_H
#define immersedBasicSymmetryFvPatchField_H

#include "immersedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class immersedBasicSymmetryFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class immersedBasicSymmetryFvPatchField
:
    public immersedValueFvPatchField<Type>
{
public:

    // Constructors

        //- Construct from patch and internal field
        immersedBasicSymmetryFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        immersedBasicSymmetryFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given
        // immersedBasicSymmetryFvPatchField<Type> onto a new patch
        immersedBasicSymmetryFvPatchField
        (
            const immersedBasicSymmetryFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        immersedBasicSymmetryFvPatchField
        (
            const immersedBasicSymmetryFvPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        immersedBasicSymmetryFvPatchField
        (
            const immersedBasicSymmetryFvPatchField<Type>&,
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
                new immersedBasicSymmetryFvPatchField<Type>(*this, iF)
            );
        }


    // Evaluation functions

        //- Update the boundary values and set the patchInternal field
        virtual void updateCoeffs();

        //- Set the internal values
        virtual void setInternal()
        {}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "immersedBasicSymmetryFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
