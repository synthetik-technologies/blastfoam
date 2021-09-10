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
    Foam::bisectionRootSolver

Description
    Find root using the bisection method

SourceFiles
    rootSolver.C

\*---------------------------------------------------------------------------*/

#ifndef bisectionRootSolver_H
#define bisectionRootSolver_H

#include "rootSolver.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class rootSolver Declaration
\*---------------------------------------------------------------------------*/

class bisectionRootSolver
:
    public rootSolver
{
public:

    //- Runtime type information
    TypeName("bisection");

    // Constructors

        //- Construct for given rootSolver
        bisectionRootSolver(const rootSystem& eqn, const dictionary& dict);

        //- Disallow default bitwise copy construction
        bisectionRootSolver(const rootSolver&) = delete;


    //- Destructor
    virtual ~bisectionRootSolver()
    {}


    // Member Functions

        //- Solve the root
        virtual scalar solve
        (
            const scalar x0,
            const label li
        ) const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const bisectionRootSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
