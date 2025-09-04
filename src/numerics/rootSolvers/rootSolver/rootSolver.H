/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    Generic class for finding the roots of a multivariate system

SourceFiles
    rootSolver.C

\*---------------------------------------------------------------------------*/

#ifndef rootSolver_H
#define rootSolver_H

#include "MultivariateEquationsFwd.H"
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

        //- Reference to scalarMultivariateEquation
        const scalarMultivariateEquation& eqns_;

        //- Saftey factor when bounds of equation are crossed
        scalar boundsFac_;

        //- Absolute convergence tolerances
        scalarField absTolerances_;

        //- Relative convergence tolerances
        scalarField relTolerances_;

        //- Used convergence tolerance
        mutable scalarField tolerances_;

        //- The maximum number of sub-steps allowed for the integration step
        label maxSteps_;

        //- Current step number
        mutable label stepi_;

        //- Current x errors
        mutable scalarField errors_;



    // Protected member function

        //- Normalise x tolerance
        void initialise(const scalarList& x) const;

        //- Is the solution converged
        bool converged
        (
            const scalarList& dx
        ) const;

        //- Is the solution converged
        bool converged
        (
            const scalarList& x0,
            const scalarList& x1
        ) const;

        //- Print information about the current step
        void printStepInformation(const scalarList& vals) const;

        //- Print information about the current step
        void printFinalInformation(const scalarList& vals) const;


public:

    //- Runtime type information
    TypeName("rootSolver");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            rootSolver,
            dictionaryUnivariate,
            (const scalarMultivariateEquation& eqns, const dictionary& dict),
            (eqns, dict)
        );
        declareRunTimeSelectionTable
        (
            autoPtr,
            rootSolver,
            dictionaryZero,
            (const scalarMultivariateEquation& eqns, const dictionary& dict),
            (eqns, dict)
        );
        declareRunTimeSelectionTable
        (
            autoPtr,
            rootSolver,
            dictionaryOne,
            (const scalarMultivariateEquation& eqns, const dictionary& dict),
            (eqns, dict)
        );


    // Constructors

        //- Construct for given rootSolver
        rootSolver
        (
            const scalarMultivariateEquation& eqns,
            const dictionary& dict
        );

        //- Copy constructor
        rootSolver
        (
            const scalarMultivariateEquation& eqns,
            const rootSolver& solver
        );

        //- Disallow default bitwise copy construction
        rootSolver(const rootSolver&) = delete;

    // Selectors

        //- Select null constructed
        static autoPtr<rootSolver> New
        (
            const scalarMultivariateEquation& eqns,
            const dictionary& dict
        );

        //- Select null constructed
        static autoPtr<rootSolver> New
        (
            const word& rootSolverType,
            const scalarMultivariateEquation& eqns,
            const dictionary& dict
        );



        //- Clone the root solver
        virtual autoPtr<rootSolver> clone
        (
            const scalarMultivariateEquation& eqns
        ) const = 0;


    //- Destructor
    virtual ~rootSolver();


    // Member Functions

        //- Return the tolerances
        inline const scalarField& tolerances() const;

        //- Return the relative tolerances
        inline const scalarField& relTolerances() const;

        //- Access the relative tolerances
        inline scalarField& relTolerances();

        //- Return the absolute tolerances
        inline const scalarField& absTolerances() const;

        //- Access the absolute tolerances
        inline scalarField& absTolerances();

        //- Return the step number
        inline label nSteps() const;

        //- Return the x errors
        inline const scalarList& errors() const;

        //- Find the roots of the equation
        virtual tmp<scalarField> solve() const;

        //- Find the roots of the equation
        virtual tmp<scalarField> solve(const scalarList& x0) const;

        //- Find the roots of the equation between xHigh and xLow
        virtual tmp<scalarField> solve
        (
            const scalarList& x0,
            const scalarList& xLow,
            const scalarList& xHigh
        ) const;

        //- Find the roots of the equation
        virtual tmp<scalarField> solve
        (
            const scalarList& x0,
            const label li
        ) const;

        //- Find the roots of the equation between xHigh and xLow
        virtual tmp<scalarField> solve
        (
            const scalarList& x0,
            const scalarList& xLow,
            const scalarList& xHigh,
            const label li
        ) const;

        //- Find the roots of the equation between xHigh and xLow
        virtual tmp<scalarField> findRoots
        (
            const scalarList& x0,
            const scalarList& xLow,
            const scalarList& xHigh,
            const label li
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rootSolverI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
