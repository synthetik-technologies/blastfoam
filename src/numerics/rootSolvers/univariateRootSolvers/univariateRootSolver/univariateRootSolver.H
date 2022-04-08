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
    Foam::univariateRootSolver

Description
    Generic class for finding a root of a univariate equation

SourceFiles
    univariateRootSolver.C

\*---------------------------------------------------------------------------*/

#ifndef univariateRootSolver_H
#define univariateRootSolver_H

#include "rootSolver.H"
#include "EquationsFwd.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class univariateRootSolver Declaration
\*---------------------------------------------------------------------------*/

class univariateRootSolver
:
    public rootSolver
{

protected:

    // Protected data

        //- Reference to scalarEquation
        const scalarEquation& eqn_;


    // Protected member functions

        //- Normalise x tolerance
        void initialise(const scalar x) const;

        //- Is the solution converged
        bool converged
        (
            const scalar dx,
            const scalar y
        ) const;

        //- Is the solution converged
        bool converged
        (
            const scalar x0,
            const scalar x1,
            const scalar y
        ) const;

        //- Print information about the current step
        void printStepInformation(const scalar val) const;

        //- Print information about the current step
        scalar printFinalInformation(const scalar val) const;


public:

    //- Runtime type information (debug defined in rootSolver)
    TypeNameNoDebug("univariateRootSolver");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            univariateRootSolver,
            dictionaryZero,
            (const scalarMultivariateEquation& eqn, const dictionary& dict),
            (eqn, dict)
        );
        declareRunTimeSelectionTable
        (
            autoPtr,
            univariateRootSolver,
            dictionaryOne,
            (const scalarMultivariateEquation& eqn, const dictionary& dict),
            (eqn, dict)
        );
        declareRunTimeSelectionTable
        (
            autoPtr,
            univariateRootSolver,
            dictionaryTwo,
            (const scalarMultivariateEquation& eqn, const dictionary& dict),
            (eqn, dict)
        );


    // Constructors

        //- Construct for given univariateRootSolver
        univariateRootSolver
        (
            const scalarMultivariateEquation& eqn,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        univariateRootSolver(const univariateRootSolver&) = delete;


    // Selectors

        //- Select from a dictionary
        static autoPtr<univariateRootSolver> New
        (
            const scalarEquation& eqn,
            const dictionary& dict
        );

        //- Select given a type name
        static autoPtr<univariateRootSolver> New
        (
            const word& rootSolverType,
            const scalarEquation& eqn,
            const dictionary& dict
        );


    //- Destructor
    virtual ~univariateRootSolver();


    // Member Functions

        //- Return the x error
        inline scalar xError() const;

        //- Return the y error
        inline scalar yError() const;

        //- Return the x tolerance
        inline scalar xTol() const;

        //- Access the x tolerance
        inline scalar& xTol();

        //- Return the absolute x tolerance
        inline scalar xAbsTol() const;

        //- Access the absolute x tolerance
        inline scalar& xAbsTol();

        //- Return the y tolerance
        inline scalar yTol() const;

        //- Access the y tolerance
        inline scalar& yTol();

        using rootSolver::solve;

        //- Find the root of the equation
        virtual scalar solve(const scalar x0) const;

        //- Find the root of the equation between xHigh and xLow
        virtual scalar solve
        (
            const scalar x0,
            const scalar xLow,
            const scalar xHigh
        ) const;

        //- Find the root of the equation
        virtual scalar solve
        (
            const scalar x0,
            const label li
        ) const;

        //- Find the root of the equation between xHigh and xLow
        virtual scalar solve
        (
            const scalar x0,
            const scalar xLow,
            const scalar xHigh,
            const label li
        ) const;

        //- Find the root of the equation between xHigh and xLow
        virtual tmp<scalarField> findRoots
        (
            const scalarField& x0,
            const scalarField& xLow,
            const scalarField& xHigh,
            const label li
        ) const;

        //- Find the root of the equation between xHigh and xLow
        virtual scalar findRoot
        (
            const scalar x0,
            const scalar xLow,
            const scalar xHigh,
            const label li
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "univariateRootSolverI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
