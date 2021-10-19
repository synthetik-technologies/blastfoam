/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    RBFFunction

Description
    Radial basis function virtual base class

Author
    Frank Bos, TU Delft.  All rights reserved.

SourceFiles
    RBFFunction.C
    newRBFFunction.C

\*---------------------------------------------------------------------------*/

#ifndef RBFFunction_H
#define RBFFunction_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "tmp.H"
#include "autoPtr.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class RBFFunction Declaration
\*---------------------------------------------------------------------------*/

class RBFFunction
{
    // Private Member Functions

        //- Disallow copy construct
        RBFFunction(const RBFFunction&);

        //- Disallow default bitwise assignment
        void operator=(const RBFFunction&);


public:

    //- Runtime type information
    TypeName("RBFFunction");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            RBFFunction,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Selectors

        //- Return a pointer to the selected RBF function
        static autoPtr<RBFFunction> New(const dictionary& dict);

        //- Return a pointer to the selected RBF function
        static autoPtr<RBFFunction> New
        (
            const word& type,
            const dictionary& dict
        );



    // Constructors

        //- Construct null
        RBFFunction()
        {}

        //- Create and return a clone
        virtual autoPtr<RBFFunction> clone() const = 0;


    // Destructor

        virtual ~RBFFunction()
        {}


    // Member Functions

        //- Return RBF weights
        virtual tmp<scalarField> weights
        (
            const vectorField& controlPoints,
            const vector& dataPoint
        ) const;

        virtual scalar evaluate(const scalar dist) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
