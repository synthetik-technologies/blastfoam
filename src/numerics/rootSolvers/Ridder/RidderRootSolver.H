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
    Foam::RidderRootSolver

Description
    Find root using the Ridder method

SourceFiles
    RidderRootSolver.C

\*---------------------------------------------------------------------------*/

#ifndef RidderRootSolver_H
#define RidderRootSolver_H

#include "rootSolver.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class RidderRootSolver Declaration
\*---------------------------------------------------------------------------*/

class RidderRootSolver
:
    public rootSolver
{
public:

    //- Runtime type information
    TypeName("Ridder");

    // Constructors

        //- Construct for given rootSolver
        RidderRootSolver(const equation& eqn, const dictionary& dict);

        //- Disallow default bitwise copy construction
        RidderRootSolver(const rootSolver&) = delete;


    //- Destructor
    virtual ~RidderRootSolver()
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
        void operator=(const RidderRootSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
