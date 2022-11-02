/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
01-11-2022 Synthetik Applied Technologies : Added Function4
-------------------------------------------------------------------------------
License
    This file is a derived work of OpenFOAM.

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
    Foam::Function4

Description
    Run-time selectable function of four variables (t, x, y, z)

SourceFiles
    Function4.C
    Function4New.C

\*---------------------------------------------------------------------------*/

#ifndef Function4_H
#define Function4_H

#include "dictionary.H"
#include "tmp.H"
#include "typeInfo.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
template<class Type> class Function4;
template<class Type> Ostream& operator<<(Ostream&, const Function4<Type>&);

/*---------------------------------------------------------------------------*\
                          Class Function4 Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Function4
:
    public tmp<Function4<Type>>::refCount
{

protected:

    // Protected data

        //- Name of entry
        const word name_;


public:

    typedef Type returnType;

    //- Runtime type information
    TypeName("Function4")

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        Function4,
        dictionary,
        (
            const word& name,
            const dictionary& dict
        ),
        (name, dict)
    );


    // Constructors

        //- Construct from name
        Function4(const word& name);

        //- Copy constructor
        Function4(const Function4<Type>& f1);

        //- Construct and return a clone
        virtual tmp<Function4<Type>> clone() const = 0;


    //- Selector
    static autoPtr<Function4<Type>> New
    (
        const word& name,
        const dictionary& dict
    );


    //- Destructor
    virtual ~Function4();


    // Member Functions

        //- Return the name of the entry
        const word& name() const;

        //- Return value as a function of two scalars
        virtual Type value
        (
            const scalar t,
            const scalar x,
            const scalar y,
            const scalar z
        ) const = 0;

        //- Return value as a function of two scalar fields
        virtual tmp<Field<Type>> value
        (
            const scalar t,
            const scalarField& x,
            const scalarField& y,
            const scalarField& z
        ) const = 0;

        //- Write in dictionary format
        virtual void write(Ostream& os) const = 0;


    // Member Operators

        //- Assignment
        void operator=(const Function4<Type>&);


    // IOstream Operators

        //- Ostream Operator
        friend Ostream& operator<< <Type>
        (
            Ostream& os,
            const Function4<Type>& func
        );
};


template<class Type>
void writeEntry(Ostream& os, const Function4<Type>& f1);


/*---------------------------------------------------------------------------*\
                       Class FieldFunction4 Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class Function4Type>
class FieldFunction4
:
    public Function4<Type>
{

public:

    // Constructors

        //- Construct from name
        FieldFunction4(const word& name);

        //- Construct and return a clone
        virtual tmp<Function4<Type>> clone() const;


    //- Destructor
    virtual ~FieldFunction4();


    // Member Functions

        //- Return value as a function of two scalars
        virtual Type value
        (
            const scalar t,
            const scalar x,
            const scalar y,
            const scalar z
        ) const = 0;

        //- Return value as a function of two scalar fields
        virtual tmp<Field<Type>> value
        (
            const scalar t,
            const scalarField& x,
            const scalarField& y,
            const scalarField& z
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeFunction4(Type)                                                    \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(Function4<Type>, 0);                   \
                                                                               \
    defineTemplateRunTimeSelectionTable(Function4<Type>, dictionary);


#define makeFunction4Type(SS, Type)                                            \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(Function4s::SS<Type>, 0);              \
                                                                               \
    Function4<Type>::adddictionaryConstructorToTable<Function4s::SS<Type>>     \
        addFunction4##SS##Type##ConstructorToTable_;


#define makeScalarFunction4(SS)                                                \
                                                                               \
    defineTypeNameAndDebug(SS, 0);                                             \
                                                                               \
    Function4<scalar>::adddictionaryConstructorToTable<SS>                     \
        addFunction4##SS##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Function4.C"
    #include "Constant4.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
