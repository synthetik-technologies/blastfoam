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
    Foam::SteffensenRootSolver

Description
    Find root using the Steffensen method

SourceFiles
    SteffensenRootSolver.C

\*---------------------------------------------------------------------------*/

#ifndef SteffensenRootSolver_H
#define SteffensenRootSolver_H

#include "rootSolver.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class SteffensenRootSolver Declaration
\*---------------------------------------------------------------------------*/

class SteffensenRootSolver
:
    public rootSolver
{
public:

    //- Runtime type information
    TypeName("Steffensen");

    // Constructors

        //- Construct for given rootSolver
        SteffensenRootSolver(const equation& eqn, const dictionary& dict);

        //- Disallow default bitwise copy construction
        SteffensenRootSolver(const rootSolver&) = delete;


    //- Destructor
    virtual ~SteffensenRootSolver()
    {}


    // Member Functions

        //- Solve the root
        virtual scalar findRoot
        (
            const scalar x0,
            const scalar xLow,
            const scalar xHigh,
            const label li
        ) const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const SteffensenRootSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
