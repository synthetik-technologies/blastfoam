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
    Foam::immersedValueFvPatchField

Description


SourceFiles
    immersedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef immersedValueFvPatchField_H
#define immersedValueFvPatchField_H

#include "volFields.H"
#include "immersedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class immersedValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class immersedValueFvPatchField
:
    public immersedFvPatchField<Type>
{
protected:

    //- Value
    Field<Type> immersedField_;


public:

    // Constructors

        //- Construct from patch and internal field
        immersedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        immersedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );

        //- Construct by mapping the given
        // immersedValueFvPatchField<Type> onto a new patch
        immersedValueFvPatchField
        (
            const immersedValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        immersedValueFvPatchField
        (
            const immersedValueFvPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        immersedValueFvPatchField
        (
            const immersedValueFvPatchField<Type>&,
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
                new immersedValueFvPatchField<Type>(*this, iF)
            );
        }

    //- Read internal value
    void readInternal(const dictionary& dict);

    // Member Functions

        // Map functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap( const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchField<Type>&, const labelList&);

        // Access

            //- Return the boundary values
            const Field<Type>& immersedField() const
            {
                return immersedField_;
            }

            //- Access the boundary values
            Field<Type>& immersedField()
            {
                return immersedField_;
            }

        // Evaluation functions

            //- Update the boundary values and set the patchInternal field
            virtual void updateCoeffs();

            //- Set the internal values
            virtual void setInternal()
            {}

            //- Add forcing to the given field
            virtual void addForcing
            (
                Field<Type>& f,
                const Field<scalar>& alphaRho,
                const Field<Type>& old,
                const Field<Type>& RHS,
                const scalar& dt
            ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "immersedValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
