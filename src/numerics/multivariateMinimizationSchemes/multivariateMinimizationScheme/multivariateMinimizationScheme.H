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
    Foam::multivariateMinimizationScheme

Description
    Abstract base-class for root system solvers

SourceFiles
    multivariateMinimizationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef multivariateMinimizationScheme_H
#define multivariateMinimizationScheme_H

#include "scalarMultivariateEquation.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class multivariateMinimizationScheme Declaration
\*---------------------------------------------------------------------------*/

class multivariateMinimizationScheme
{

protected:

    // Protected data

        //- Reference to scalarEquation
        const scalarMultivariateEquation& eqns_;

        //- Convergence tolerance
        scalarField xTolerances_;
        scalarField yTolerances_;

        //- The maximum number of sub-steps allowed for the integration step
        label maxSteps_;

        //- Current step number
        mutable label stepi_;

        //- Current errors
        mutable scalarField xErrors_;
        mutable scalarField yErrors_;

        //- Line search reduction coefficient
        scalar tau_;

        //- Line search
        scalar lineSearch
        (
            const scalarList& x0,
            const scalarList& grad,
            const label li,
            scalarList& fx
        ) const;

        //- Return the L2 norm the given list
        scalar norm(const scalarList& lst) const;

        //- Is the solution converged
        bool converged
        (
            const scalarList& xError, 
            const scalarList& yError
        ) const;

        //- Print information about the current step
        void printStepInformation(const scalarList& vals) const;

        //- Print information about the current step
        void printFinalInformation() const;


public:

    //- Runtime type information
    TypeName("multivariateMinimizationScheme");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            multivariateMinimizationScheme,
            dictionaryZero,
            (const scalarMultivariateEquation& eqns, const dictionary& dict),
            (eqns, dict)
        );
        declareRunTimeSelectionTable
        (
            autoPtr,
            multivariateMinimizationScheme,
            dictionaryOne,
            (const scalarMultivariateEquation& eqns, const dictionary& dict),
            (eqns, dict)
        );
        declareRunTimeSelectionTable
        (
            autoPtr,
            multivariateMinimizationScheme,
            dictionaryTwo,
            (const scalarMultivariateEquation& eqns, const dictionary& dict),
            (eqns, dict)
        );


    // Constructors

        //- Construct for given multivariateMinimizationScheme
        multivariateMinimizationScheme
        (
            const scalarMultivariateEquation& eqns,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        multivariateMinimizationScheme(const multivariateMinimizationScheme&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<multivariateMinimizationScheme> New
        (
            const scalarMultivariateEquation& eqn,
            const dictionary& dict
        );


    //- Destructor
    virtual ~multivariateMinimizationScheme()
    {}


    // Member Functions

        //- Return access to the tolerance
        inline scalarList& xTolerances();
        inline scalarList& yTolerances();

        //- Return the step number
        inline label nSteps() const;

        //- Return the errors
        inline scalarList& xErrors() const;
        inline scalarList& yErrors() const;

        //- Find the root of the equation
        virtual tmp<scalarField> solve() const;

        //- Find the root of the equation
        virtual tmp<scalarField> solve(const scalarList& x0) const;

        //- Find the root of the equation between xHigh and xLow
        virtual tmp<scalarField> solve
        (
            const scalarList& x0,
            const scalarList& xLow,
            const scalarList& xHigh
        ) const;

        //- Find the root of the equation
        virtual tmp<scalarField> solve
        (
            const scalarList& x0,
            const label li
        ) const;

        //- Find the root of the equation between xHigh and xLow
        virtual tmp<scalarField> solve
        (
            const scalarList& x0,
            const scalarList& xLow,
            const scalarList& xHigh,
            const label li
        ) const;

        //- Find the root of the equation between xHigh and xLow
        virtual tmp<scalarField> minimize
        (
            const scalarList& x0,
            const scalarList& xLow,
            const scalarList& xHigh,
            const label li
        ) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multivariateMinimizationScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multivariateMinimizationSchemeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
