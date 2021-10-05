/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::bisectionMinimizationScheme

Description
    Find root using the bisection method

SourceFiles
    bisectionMinimizationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef bisectionMinimizationScheme_H
#define bisectionMinimizationScheme_H

#include "minimizationScheme.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class bisectionMinimizationScheme Declaration
\*---------------------------------------------------------------------------*/

class bisectionMinimizationScheme
:
    public minimizationScheme
{
public:

    //- Runtime type information
    TypeName("bisection");

    // Constructors

        //- Construct for given minimizationScheme
        bisectionMinimizationScheme(const scalarEquation& eqn, const dictionary& dict);

        //- Disallow default bitwise copy construction
        bisectionMinimizationScheme(const minimizationScheme&) = delete;


    //- Destructor
    virtual ~bisectionMinimizationScheme()
    {}


    // Member Functions

        //- Solve the root
        using minimizationScheme::solve;
        virtual scalar solve
        (
            const scalar x0,
            const scalar x1,
            const scalar x2,
            const label li
        ) const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const bisectionMinimizationScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
