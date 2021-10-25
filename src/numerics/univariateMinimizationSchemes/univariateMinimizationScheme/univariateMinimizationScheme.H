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
    Foam::univariateMinimizationScheme

Description
    Abstract base-class for root system solvers

SourceFiles
    univariateMinimizationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef univariateMinimizationScheme_H
#define univariateMinimizationScheme_H

#include "equation.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class univariateMinimizationScheme Declaration
\*---------------------------------------------------------------------------*/

class univariateMinimizationScheme
{

protected:

    // Protected data

        //- Reference to scalarEquation
        const equation& eqn_;

        //- Convergence tolerance
        scalar tolerance_;

        //- The maximum number of sub-steps allowed for the integration step
        label maxSteps_;

        //- Number of location to pre sample
        label nSamples_;

        //- Current step number
        mutable label stepi_;

        //- Current error
        mutable scalar error_;

        //- Is the solution converged
        bool converged(const scalar error) const;

        //- Print information about the current step
        void printStepInformation(const scalar val) const;

        //- Print information about the current step
        scalar printFinalInformation(const scalar val) const;

        //- Sample and return the interval with the smallest value
        void sample(scalar&, scalar&, const label) const;

public:

    //- Runtime type information
    TypeName("univariateMinimizationScheme");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            univariateMinimizationScheme,
            dictionaryZero,
            (const equation& eqn, const dictionary& dict),
            (eqn, dict)
        );
        declareRunTimeSelectionTable
        (
            autoPtr,
            univariateMinimizationScheme,
            dictionaryOne,
            (const equation& eqn, const dictionary& dict),
            (eqn, dict)
        );
        declareRunTimeSelectionTable
        (
            autoPtr,
            univariateMinimizationScheme,
            dictionaryTwo,
            (const equation& eqn, const dictionary& dict),
            (eqn, dict)
        );


    // Constructors

        //- Construct for given univariateMinimizationScheme
        univariateMinimizationScheme(const equation& eqn, const dictionary& dict);

        //- Disallow default bitwise copy construction
        univariateMinimizationScheme(const univariateMinimizationScheme&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<univariateMinimizationScheme> New
        (
            const equation& eqn,
            const dictionary& dict
        );


    //- Destructor
    virtual ~univariateMinimizationScheme()
    {}


    // Member Functions

        //- Return access to the tolerance
        inline scalar& tolerance();

        //- Return the step number
        inline label nSteps() const;

        //- Return the error
        inline scalar error() const;

        //- Find the root of the equation
        virtual scalar solve() const;

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
        virtual scalar minimize
        (
            const scalar x0,
            const scalar xLow,
            const scalar xHigh,
            const label li
        ) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const univariateMinimizationScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "univariateMinimizationSchemeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
