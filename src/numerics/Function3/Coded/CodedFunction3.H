/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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
    Foam::Function3s::Coded

Description
    Constructs a dynamically compiled function of three variables.
    Helper functions have been added to easily create helper function1
    or function2 classes using Func1(Name, var) and Func2(Name, var).
    the macros will return a non-constant reference to a runTime
    selectable function of the type specified in the Name dictionary.
    Variables in the code dictionary are updated when the object is made.

See also
    Foam::dynamicCode
    Foam::functionEntries::codeStream

SourceFiles
    CodedFunction3I.H
    CodedFunction3.C

\*---------------------------------------------------------------------------*/

#ifndef CodedFunction3_H
#define CodedFunction3_H

#include "Function3.H"
#include "codedBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function3s
{

/*---------------------------------------------------------------------------*\
                            Class coded Declaration
\*---------------------------------------------------------------------------*/

class coded
{
public:

    ClassNameNoDebug("Function3");
};


/*---------------------------------------------------------------------------*\
                            Class Coded Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Coded
:
    public Function3<Type>,
    public codedBase
{
    // Private Data

        //- Saved dictionary
        const dictionary& dict_;

        //- The dynamically generated Function3 pointer
        mutable autoPtr<Function3<Type>> redirectFunction3Ptr_;


    // Private Member Functions

        //- Get the keywords associated with source code
        virtual wordList codeKeys() const;

        //- Adapt the context for the current object
        virtual void prepare(dynamicCode&, const dynamicCodeContext&) const;

        //- Clear the ptr to the redirected object
        virtual void clearRedirect() const;

        //- Compile, link and return the now coded Function3
        autoPtr<Function3<Type>> compileNew();

        //- Expand strings in the code entry
        const dictionary& expandCodeDict(const dictionary& dict) const;


public:

    // Runtime type information
    TypeName("coded");


    // Constructors

        //- Construct from name and dictionary
        Coded
        (
            const word& name,
            const dictionary& dict
        );

        //- Copy constructor
        Coded(const Coded<Type>& cf1);

        //- Construct and return a clone
        virtual tmp<Function3<Type>> clone() const;


    //- Destructor
    virtual ~Coded();


    // Member Functions

        //- Return value as a function of two scalars
        virtual inline Type value
        (
            const scalar x,
            const scalar y,
            const scalar z
        ) const;

        //- Return value as a function of two scalar fields
        virtual tmp<Field<Type>> value
        (
            const scalarField& x,
            const scalarField& y,
            const scalarField& z
        ) const;

        //- Write data to dictionary stream
        virtual void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Coded<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function3s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CodedFunction3I.H"

#ifdef NoRepository
    #include "CodedFunction3.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
