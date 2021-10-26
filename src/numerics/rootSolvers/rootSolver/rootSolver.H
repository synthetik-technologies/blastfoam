/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
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

#include "MultivariateEquation.H"
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
        const multivariateEquation<scalar>& eqns_;

        //- Convergence tolerance
        scalarField tolerances_;

        //- The maximum number of sub-steps allowed for the integration step
        label maxSteps_;

        //- Current step number
        mutable label stepi_;

        //- Current errors
        mutable scalarField errors_;

        //- Is the solution converged
        bool converged(const scalarField& errors) const;

        //- Print information about the current step
        void printStepInformation(const scalarField& vals) const;

        //- Print information about the current step
        void printFinalInformation() const;

public:

    //- Runtime type information
    TypeName("rootSolver");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            rootSolver,
            dictionaryUnivariate,
            (const multivariateEquation<scalar>& eqns, const dictionary& dict),
            (eqns, dict)
        );
        declareRunTimeSelectionTable
        (
            autoPtr,
            rootSolver,
            dictionaryZero,
            (const multivariateEquation<scalar>& eqns, const dictionary& dict),
            (eqns, dict)
        );
        declareRunTimeSelectionTable
        (
            autoPtr,
            rootSolver,
            dictionaryOne,
            (const multivariateEquation<scalar>& eqns, const dictionary& dict),
            (eqns, dict)
        );


    // Constructors

        //- Construct for given rootSolver
        rootSolver
        (
            const multivariateEquation<scalar>& eqns,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        rootSolver(const rootSolver&) = delete;

    // Selectors

        //- Select null constructed
        static autoPtr<rootSolver> New
        (
            const multivariateEquation<scalar>& eqns,
            const dictionary& dict
        );


    //- Destructor
    virtual ~rootSolver();


    // Member Functions

        //- Return access to the tolerance
        inline const scalarField& tolerances() const;

        //- Return the step number
        inline label nSteps() const;

        //- Return the errors
        inline const scalarList& errors() const;

        //- Find the roots of the equation
        virtual tmp<scalarField> solve() const;

        //- Find the roots of the equation
        virtual tmp<scalarField> solve(const scalarField& x0) const;

        //- Find the roots of the equation between xHigh and xLow
        virtual tmp<scalarField> solve
        (
            const scalarField& x0,
            const scalarField& xLow,
            const scalarField& xHigh
        ) const;

        //- Find the roots of the equation
        virtual tmp<scalarField> solve
        (
            const scalarField& x0,
            const label li
        ) const;

        //- Find the roots of the equation between xHigh and xLow
        virtual tmp<scalarField> solve
        (
            const scalarField& x0,
            const scalarField& xLow,
            const scalarField& xHigh,
            const label li
        ) const;

        //- Find the roots of the equation between xHigh and xLow
        virtual tmp<scalarField> findRoots
        (
            const scalarField& x0,
            const scalarField& xLow,
            const scalarField& xHigh,
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
