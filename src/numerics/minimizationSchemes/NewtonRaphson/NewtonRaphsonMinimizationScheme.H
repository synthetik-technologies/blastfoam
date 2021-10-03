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
    Foam::NewtonRaphsonMinimizationScheme

Description
    Find root using the NewtonRaphson method

SourceFiles
    NewtonRaphsonMinimizationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef NewtonRaphsonMinimizationScheme_H
#define NewtonRaphsonMinimizationScheme_H

#include "minimizationScheme.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class NewtonRaphsonMinimizationScheme Declaration
\*---------------------------------------------------------------------------*/

class NewtonRaphsonMinimizationScheme
:
    public minimizationScheme
{
public:

    //- Runtime type information
    TypeName("NewtonRaphson");

    // Constructors

        //- Construct for given minimizationScheme
        NewtonRaphsonMinimizationScheme(const scalarEquation& eqn, const dictionary& dict);

        //- Disallow default bitwise copy construction
        NewtonRaphsonMinimizationScheme(const minimizationScheme&) = delete;


    //- Destructor
    virtual ~NewtonRaphsonMinimizationScheme()
    {}


    // Member Functions

        //- Solve the root
        virtual scalar solve
        (
            const scalar x0,
            const scalar x1,
            const scalar x2,
            const label li
        ) const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const NewtonRaphsonMinimizationScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
