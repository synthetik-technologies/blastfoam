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
    Foam::Function4s::Scale

Description
    Function4 which scales a given 'value' function by a 'scale' scalar function
    and scales the 'x','y', and 'z' arguments of the 'value' and 'scale' functions
    by the optional 'xScale', 'yScale', and 'zScale' scalar functions.

SourceFiles
    Scale.C

\*---------------------------------------------------------------------------*/

#ifndef Scale4_H
#define Scale4_H

#include "Function4.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function4s
{

/*---------------------------------------------------------------------------*\
                            Class Scale Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Scale
:
    public FieldFunction4<Type, Scale<Type>>
{
    // Private Data

        //- Scalar scaling function
        autoPtr<Function4<scalar>> scale_;

        //- Time scaling function
        autoPtr<Function1<scalar>> tScale_;

        //- Argument scaling function
        autoPtr<Function1<scalar>> xScale_;

        //- Argument scaling function
        autoPtr<Function1<scalar>> yScale_;

        //- Argument scaling function
        autoPtr<Function1<scalar>> zScale_;

        //- Value function
        autoPtr<Function4<Type>> value_;


    // Private Member Functions

        //- Read the coefficients from the given dictionary
        void read
        (
            const dictionary& dict,
            const unitConversions& units
        );


public:

    // Runtime type information
    TypeName("scale");


    // Constructors

        //- Construct from name and dictionary
        Scale
        (
            const word& name,
            const unitConversions& units,
            const dictionary& dict
        );

        //- Copy constructor
        Scale(const Scale<Type>& se);


    //- Destructor
    virtual ~Scale();


    // Member Functions

        //- Return value
        virtual inline Type value
        (
            const scalar t,
            const scalar x,
            const scalar y,
            const scalar z
        ) const;

        //- Write data to dictionary stream
        virtual void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Scale<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function4s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Scale4I.H"

#ifdef NoRepository
    #include "Scale4.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
