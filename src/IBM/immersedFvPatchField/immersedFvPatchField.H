/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::immersedFvPatchField

Description
    This boundary condition provides an 'immersed' condition for reduced
    dimensions cases, i.e. 1- and 2-D geometries.  Apply this condition to
    patches whose normal is aligned to geometric directions that do not
    constitue solution directions.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            immersed;
    }
    \endverbatim

SourceFiles
    immersedFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef immersedFvPatchField_H
#define immersedFvPatchField_H

#include "fvPatchField.H"
#include "immersedFvPatch.H"
#include "immersedBoundaryObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class immersedFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class immersedFvPatchField
:
    public fvPatchField<Type>
{
protected:

    // Reference to the underlying immersed object
    immersedBoundaryObject& object_;


public:

    //- Runtime type information
    TypeName(immersedFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        immersedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        immersedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given immersedFvPatchField onto a new patch
        immersedFvPatchField
        (
            const immersedFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        immersedFvPatchField(const immersedFvPatchField<Type>&) = delete;

        //- Copy constructor setting internal field reference
        immersedFvPatchField
        (
            const immersedFvPatchField<Type>&,
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
                new immersedFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap( const fvPatchFieldMapper&)
            {}

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchField<Type>&, const labelList&)
            {}


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            //  Sets Updated to true
            virtual void updateCoeffs();

            //- Return internal field next to patch as patch field
//             virtual tmp<Field<Type>> patchInternalField() const;

            //- Return internal field next to patch as patch field
//             virtual void patchInternalField(Field<Type>&) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const
            {
                return tmp<Field<Type>>(new Field<Type>(0));
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const
            {
                return tmp<Field<Type>>(new Field<Type>(0));
            }

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            tmp<Field<Type>> gradientInternalCoeffs() const
            {
                return tmp<Field<Type>>(new Field<Type>(0));
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            tmp<Field<Type>> gradientBoundaryCoeffs() const
            {
                return tmp<Field<Type>>(new Field<Type>(0));
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "immersedFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
