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
    Foam::Function4s::Constant

Description
    Templated function of three variables that returns a constant value.

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

#ifndef Constant4_H
#define Constant4_H

#include "Function4.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function4s
{

/*---------------------------------------------------------------------------*\
                          Class Constant Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Constant
:
    public FieldFunction4<Type, Constant<Type>>
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
        //  Reads the constant value without the Function4 type
        //  for backward compatibility
        Constant(const word& name, Istream& is);

        //- Copy constructor
        Constant(const Constant<Type>& cnst);

        //- Construct and return a clone
        virtual tmp<Function4<Type>> clone() const
        {
            return tmp<Function4<Type>>(new Constant<Type>(*this));
        }


    //- Destructor
    virtual ~Constant();


    // Member Functions

        //- Return constant
        virtual inline Type value
        (
            const scalar t,
            const scalar x,
            const scalar y,
            const scalar z
        ) const;

        //- Return constant field
        virtual inline tmp<Field<Type>> value
        (
            const scalar t,
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

} // End namespace Function4s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Constant4I.H"

#ifdef NoRepository
    #include "Constant4.C"
    #include "Function4New.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
