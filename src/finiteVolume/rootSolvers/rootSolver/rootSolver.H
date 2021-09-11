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
    Foam::rootSolver

Description
    Abstract base-class for root system solvers

SourceFiles
    rootSolver.C

\*---------------------------------------------------------------------------*/

#ifndef rootSolver_H
#define rootSolver_H

#include "scalarEquation.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class rootSolver Declaration
\*---------------------------------------------------------------------------*/

class rootSolver
{

protected:

    // Protected data

        //- Reference to scalarEquation
        const scalarEquation& eqn_;

        //- Convergence tolerance
        scalar tolerance_;

        //- The maximum number of sub-steps allowed for the integration step
        label maxSteps_;

        //- Current step number
        mutable label stepi_;

        //- Current error
        mutable scalar error_;

        //- Is the solution converged
        bool converged(const scalar error) const;

public:

    //- Runtime type information
    TypeName("rootSolver");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            rootSolver,
            dictionary,
            (const scalarEquation& eqn, const dictionary& dict),
            (eqn, dict)
        );


    // Constructors

        //- Construct for given rootSolver
        rootSolver(const scalarEquation& eqn, const dictionary& dict);

        //- Disallow default bitwise copy construction
        rootSolver(const rootSolver&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<rootSolver> New
        (
            const scalarEquation& eqn,
            const dictionary& dict
        );


    //- Destructor
    virtual ~rootSolver()
    {}


    // Member Functions

        //- Return access to the tolerance
        inline scalar& tolerance();

        //- Return the step number
        inline label nSteps() const;

        //- Return the error
        inline scalar error() const;

        //- Find the root of the equation
        //- Solve the root
        virtual scalar solve
        (
            const scalar x0,
            const label li
        ) const = 0;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const rootSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rootSolverI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
