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
    Foam::Function3s::Component

Description
    Function3 which uses a linear combination of f(x) + g(y) + h(z)

SourceFiles
    Component.C

\*---------------------------------------------------------------------------*/

#ifndef Component3_H
#define Component3_H

#include "Function3.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function3s
{

/*---------------------------------------------------------------------------*\
                            Class Component Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Component
:
    public FieldFunction3<Type, Component<Type>>
{
    // Private Data

        //- Argument scaling function
        autoPtr<Function1<Type>> xFunc_;

        //- Argument scaling function
        autoPtr<Function1<Type>> yFunc_;

        //- Argument scaling function
        autoPtr<Function1<Type>> zFunc_;



    // Private Member Functions

        //- Read the coefficients from the given dictionary
        void read
        (
            const dictionary& dict,
            const unitConversions& units
        );


public:

    // Runtime type information
    TypeName("component");


    // Constructors

        //- Construct from name and dictionary
        Component
        (
            const word& name,
            const unitConversions& units,
            const dictionary& dict
        );

        //- Copy constructor
        Component(const Component<Type>& se);


    //- Destructor
    virtual ~Component();


    // Member Functions

        //- Return value
        virtual inline Type value
        (
            const scalar x,
            const scalar y,
            const scalar z
        ) const;

        //- Write data to dictionary stream
        virtual void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Component<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function3s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Component3I.H"

#ifdef NoRepository
    #include "Component3.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
