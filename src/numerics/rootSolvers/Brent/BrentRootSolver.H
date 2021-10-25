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
    Foam::BrentRootSolver

Description
    Find root using the Brent method

SourceFiles
    BrentRootSolver.C

\*---------------------------------------------------------------------------*/

#ifndef BrentRootSolver_H
#define BrentRootSolver_H

#include "rootSolver.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class BrentRootSolver Declaration
\*---------------------------------------------------------------------------*/

class BrentRootSolver
:
    public rootSolver
{
public:

    //- Runtime type information
    TypeName("Brent");

    // Constructors

        //- Construct for given rootSolver
        BrentRootSolver(const equation& eqn, const dictionary& dict);

        //- Disallow default bitwise copy construction
        BrentRootSolver(const rootSolver&) = delete;


    //- Destructor
    virtual ~BrentRootSolver()
    {}


    // Member Functions

        //- Solve the root
        virtual scalar findRoot
        (
            const scalar x0,
            const scalar x1,
            const scalar x2,
            const label li
        ) const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const BrentRootSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
