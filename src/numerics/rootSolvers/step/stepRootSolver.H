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
    Foam::stepRootSolver

Description
    Find root using the step method

SourceFiles
    stepRootSolver.C

\*---------------------------------------------------------------------------*/

#ifndef stepRootSolver_H
#define stepRootSolver_H

#include "rootSolver.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class stepRootSolver Declaration
\*---------------------------------------------------------------------------*/

class stepRootSolver
:
    public rootSolver
{
// Private data

    //- Initial step size
    scalar dx_;


public:

    //- Runtime type information
    TypeName("step");

    // Constructors

        //- Construct for given rootSolver
        stepRootSolver(const equation& eqn, const dictionary& dict);

        //- Construct for given rootSolver, lower bound, and step size
        stepRootSolver
        (
            const equation& eqn,
            const scalar dx
        );


        //- Disallow default bitwise copy construction
        stepRootSolver(const rootSolver&) = delete;


    //- Destructor
    virtual ~stepRootSolver()
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
        void operator=(const stepRootSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
