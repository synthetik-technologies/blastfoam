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
    Foam::Function3s::Constant

Description
    Templated function of two variables that returns a constant value.

    Usage - for entry \<name\> returning the value <value>:
    \verbatim
        <name>     constant <value>;
    \endverbatim

    or simply

    \verbatim
        <name>     <value>;
    \endverbatim

SourceFiles
    Constant.C

\*---------------------------------------------------------------------------*/

#ifndef Constant3_H
#define Constant3_H

#include "Function3.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function3s
{

/*---------------------------------------------------------------------------*\
                          Class Constant Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Constant
:
    public FieldFunction3<Type, Constant<Type>>
{
    // Private Data

        //- Constant value
        Type value_;


public:

    // Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from name and value
        Constant(const word& name, const Type& val);

        //- Construct from name and dictionary
        Constant(const word& name, const dictionary& dict);

        //- Construct from name and Istream
        //  Reads the constant value without the Function3 type
        //  for backward compatibility
        Constant(const word& name, Istream& is);

        //- Copy constructor
        Constant(const Constant<Type>& cnst);

        //- Construct and return a clone
        virtual tmp<Function3<Type>> clone() const
        {
            return tmp<Function3<Type>>(new Constant<Type>(*this));
        }


    //- Destructor
    virtual ~Constant();


    // Member Functions

        //- Return constant
        virtual inline Type value
        (
            const scalar x,
            const scalar y,
            const scalar z
        ) const;

        //- Return constant field
        virtual inline tmp<Field<Type>> value
        (
            const scalarField& x,
            const scalarField& y,
            const scalarField& z
        ) const;

        //- Write in dictionary format
        virtual void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Constant<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function3s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Constant3I.H"

#ifdef NoRepository
    #include "Constant3.C"
    #include "Function3New.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
