/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::valueFePatchField

Description
    Foam::valueFePatchField

SourceFiles
    valueFePatchField.C

\*---------------------------------------------------------------------------*/

#ifndef valueFePatchField_H
#define valueFePatchField_H

#include "fePatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class valueFePatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class valueFePatchField
:
    public fePatchField<Type>,
    public Field<Type>
{
    // Private Member Functions

        void checkFieldSize() const;


public:

    //- Runtime type information
    TypeName("value");


    // Constructors

        //- Construct from patch and internal field
        valueFePatchField
        (
            const fePatch&,
            const DimensionedField<Type, feMesh>&
        );

        //- Construct from patch, internal field and dictionary
        valueFePatchField
        (
            const fePatch&,
            const DimensionedField<Type, feMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        valueFePatchField
        (
            const valueFePatchField<Type>&,
            const fePatch&,
            const DimensionedField<Type, feMesh>&,
            const fePatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        valueFePatchField
        (
            const valueFePatchField<Type>&,
            const DimensionedField<Type, feMesh>&
        );

        //- Disallow copy without setting internal field reference
        valueFePatchField(const valueFePatchField<Type>&) = delete;

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<fePatchField<Type>> clone
        (
            const DimensionedField<Type, feMesh>& iF
        ) const
        {
            return autoPtr<fePatchField<Type>>
            (
                new valueFePatchField<Type>
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Access

            //- Return size
            label size() const
            {
                return Field<Type>::size();
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fePatchFieldMapper&);

            //- Reverse map the given fePatchField onto this fePatchField
            //  Used to reconstruct fields
            virtual void rmap(const fePatchField<Type>&, const labelList&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );


        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        // Assignment operators

            virtual void operator=
            (
                const valueFePatchField<Type>&
            );

            virtual void operator=
            (
                const fePatchField<Type>&
            );

            virtual void operator=(const Field<Type>&);
            virtual void operator=(const Type&);


        // Force an assignment irrespective of form of patch

            virtual void operator==
            (
                const valueFePatchField<Type>&
            );

            virtual void operator==
            (
                const fePatchField<Type>&
            );

            virtual void operator==(const Field<Type>&);
            virtual void operator==(const Type&);
};


// This function is added to override the hack in fePatchField.H
// which enables simple backward compatibility with versions using
// referenceLevel in GeometicField
template<class Type>
tmp<Field<Type>> operator+
(
    const valueFePatchField<Type>& vsppf,
    const Type& t
)
{
    return static_cast<const Field<Type>&>(vsppf) + t;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "valueFePatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
