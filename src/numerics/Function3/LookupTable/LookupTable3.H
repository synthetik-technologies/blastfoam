/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
03-12-2021 Synthetik Applied Technologies : Added Function3
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
    Foam::Function3s::LookupTable

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
    LookupTable3.C

\*---------------------------------------------------------------------------*/

#ifndef LookupTable3_H
#define LookupTable3_H

#include "Function3.H"
#include "lookupTable3D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function3s
{

/*---------------------------------------------------------------------------*\
                          Class LookupTable Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class LookupTable
:
    public FieldFunction3<Type, LookupTable<Type>>
{
    // Private Data

        //- LookupTable value
        lookupTable3D<Type> table_;

        //- Units
        unitConversions units_;

public:

    // Runtime type information
    TypeName("lookupTable");


    // Constructors

        //- Construct from name and dictionary
        LookupTable
        (
            const word& name,
            const unitConversions& units,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual tmp<Function3<Type>> clone() const
        {
            return tmp<Function3<Type>>(new LookupTable<Type>(*this));
        }


    //- Destructor
    virtual ~LookupTable();


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
        void operator=(const LookupTable<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function3s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LookupTable3I.H"

#ifdef NoRepository
    #include "LookupTable3.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
