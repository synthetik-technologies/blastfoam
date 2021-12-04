/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::Function3

Description
    Run-time selectable function of two variables

    with many options provided from simple constant values to complex
    functions, interpolated tabulated data etc. etc.

SourceFiles
    Function3.C
    Function3New.C

\*---------------------------------------------------------------------------*/

#ifndef Function3_H
#define Function3_H

#include "dictionary.H"
#include "tmp.H"
#include "typeInfo.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
template<class Type> class Function3;
template<class Type> Ostream& operator<<(Ostream&, const Function3<Type>&);

/*---------------------------------------------------------------------------*\
                          Class Function3 Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Function3
:
    public tmp<Function3<Type>>::refCount
{

protected:

    // Protected data

        //- Name of entry
        const word name_;


public:

    typedef Type returnType;

    //- Runtime type information
    TypeName("Function3")

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        Function3,
        dictionary,
        (
            const word& name,
            const dictionary& dict
        ),
        (name, dict)
    );


    // Constructors

        //- Construct from name
        Function3(const word& name);

        //- Copy constructor
        Function3(const Function3<Type>& f1);

        //- Construct and return a clone
        virtual tmp<Function3<Type>> clone() const = 0;


    //- Selector
    static autoPtr<Function3<Type>> New
    (
        const word& name,
        const dictionary& dict
    );


    //- Destructor
    virtual ~Function3();


    // Member Functions

        //- Return the name of the entry
        const word& name() const;

        //- Return value as a function of two scalars
        virtual Type value
        (
            const scalar x,
            const scalar y,
            const scalar z
        ) const = 0;

        //- Return value as a function of two scalar fields
        virtual tmp<Field<Type>> value
        (
            const scalarField& x,
            const scalarField& y,
            const scalarField& z
        ) const = 0;

        //- Write in dictionary format
        virtual void write(Ostream& os) const = 0;


    // Member Operators

        //- Assignment
        void operator=(const Function3<Type>&);


    // IOstream Operators

        //- Ostream Operator
        friend Ostream& operator<< <Type>
        (
            Ostream& os,
            const Function3<Type>& func
        );
};


template<class Type>
void writeEntry(Ostream& os, const Function3<Type>& f1);


/*---------------------------------------------------------------------------*\
                       Class FieldFunction3 Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class Function3Type>
class FieldFunction3
:
    public Function3<Type>
{

public:

    // Constructors

        //- Construct from name
        FieldFunction3(const word& name);

        //- Construct and return a clone
        virtual tmp<Function3<Type>> clone() const;


    //- Destructor
    virtual ~FieldFunction3();


    // Member Functions

        //- Return value as a function of two scalars
        virtual Type value
        (
            const scalar x,
            const scalar y,
            const scalar z
        ) const = 0;

        //- Return value as a function of two scalar fields
        virtual tmp<Field<Type>> value
        (
            const scalarField& x,
            const scalarField& y,
            const scalarField& z
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeFunction3(Type)                                                    \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(Function3<Type>, 0);                   \
                                                                               \
    defineTemplateRunTimeSelectionTable(Function3<Type>, dictionary);


#define makeFunction3Type(SS, Type)                                            \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(Function3s::SS<Type>, 0);              \
                                                                               \
    Function3<Type>::adddictionaryConstructorToTable<Function3s::SS<Type>>     \
        addFunction3##SS##Type##ConstructorToTable_;


#define makeScalarFunction3(SS)                                                \
                                                                               \
    defineTypeNameAndDebug(SS, 0);                                             \
                                                                               \
    Function3<scalar>::adddictionaryConstructorToTable<SS>                     \
        addFunction3##SS##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Function3.C"
    #include "Constant3.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
