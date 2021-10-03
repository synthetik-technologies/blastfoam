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
    Foam::multivariateRootSolver

Description
    Abstract base-class for multivariateRoot system solvers

SourceFiles
    multivariateRootSolver.C

\*---------------------------------------------------------------------------*/

#ifndef multivariateRootSolver_H
#define multivariateRootSolver_H

#include "scalarMultivariateEquation.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class multivariateRootSolver Declaration
\*---------------------------------------------------------------------------*/

class multivariateRootSolver
{

protected:

    // Protected data

        //- Reference to scalarMultivariateEquation
        const scalarMultivariateEquation& eqns_;

        //- Size of the ODESystem (adjustable)
        mutable label n_;

        //- Convergence tolerance
        scalar tolerance_;

        //- The maximum number of sub-steps allowed for the integration step
        label maxSteps_;

        //- Current step number
        mutable label stepi_;

        //- Current errors
        mutable scalarField errors_;

        //- Current error
        mutable scalar error_;

        //- Is the solution converged
        bool converged(const scalarField& errors) const;

public:

    //- Runtime type information
    TypeName("multivariateRootSolver");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            multivariateRootSolver,
            dictionaryZero,
            (const scalarMultivariateEquation& eqns, const dictionary& dict),
            (eqns, dict)
        );
        declareRunTimeSelectionTable
        (
            autoPtr,
            multivariateRootSolver,
            dictionaryOne,
            (const scalarMultivariateEquation& eqns, const dictionary& dict),
            (eqns, dict)
        );


    // Constructors

        //- Construct for given multivariateRootSolver
        multivariateRootSolver
        (
            const scalarMultivariateEquation& eqns,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        multivariateRootSolver(const multivariateRootSolver&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<multivariateRootSolver> New
        (
            const scalarMultivariateEquation& eqns,
            const dictionary& dict
        );


    //- Destructor
    virtual ~multivariateRootSolver()
    {}


    // Member Functions

        //- Return access to the tolerance
        inline scalar& tolerance();

        //- Return the step number
        inline label nSteps() const;

        //- Return the errors
        inline const scalarField& errors() const;

        //- Return the error
        inline scalar error() const;

        //- Find the multivariateRoot of the equation
        //- Solve the multivariateRoot
        virtual tmp<scalarField> solve
        (
            const scalarField& x0,
            const label li
        ) const = 0;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multivariateRootSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multivariateRootSolverI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
