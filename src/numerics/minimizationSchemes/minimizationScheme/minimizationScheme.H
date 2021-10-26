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
    Foam::minimizationScheme

Description
    General class for minimizing multivariate, scalar equations

SourceFiles
    minimizationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef minimizationScheme_H
#define minimizationScheme_H

#include "ScalarEquation.H"
#include "Switch.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class minimizationScheme Declaration
\*---------------------------------------------------------------------------*/

class minimizationScheme
{

protected:

    // Protected data

        //- Reference to scalarEquation
        const scalarEquation& eqns_;

        //- Convergence tolerance
        scalarField tolerances_;

        //- The maximum number of sub-steps allowed for the integration step
        label maxSteps_;

        //- Current step number
        mutable label stepi_;

        //- Current errors
        mutable scalarField errors_;

        //- Number of pre-samples in each direction
        labelList nSamples_;

        //- Is the error normalized
        Switch normalize_;

        //- Line search reduction coefficient
        scalar tau_;

        //- Initial line search parameter
        virtual scalar alpha
        (
            const scalarField& grad,
            const scalarField& gradOld
        ) const;

        //- Line search
        scalar lineSearch
        (
            const scalarField& x0,
            const scalarField& grad,
            const scalarField& gradOld,
            const label li,
            scalar& fx
        ) const;

        //- Line search
        scalar lineSearch
        (
            const scalarField& x0,
            const scalarField& grad,
            const label li,
            scalar& fx
        ) const;

        //- Return the L2 norm the given list
        scalar norm(const scalarList& lst) const;

        //- Return lst^T * lst
        scalar inner(const scalarList& lst1, const scalarList& lst2) const;

        //- Is the solution converged
        bool converged(const scalarList& error) const;

        //- Print information about the current step
        void printStepInformation(const scalarList& vals) const;

        //- Print information about the current step
        void printFinalInformation() const;

        //- Sample and return the interval with the smallest value
        void sample
        (
            scalarField&,
            scalarField&,
            labelList& xBest,
            scalar& yBest,
            labelList& indicies,
            const label,
            const label diri
        ) const;

        //- Name of the error
        virtual word errorName() const
        {
            return "error";
        }

public:

    //- Runtime type information
    TypeName("minimizationScheme");

    // Declare run-time constructor selection tables
        declareRunTimeSelectionTable
        (
            autoPtr,
            minimizationScheme,
            dictionaryUnivariate,
            (const scalarEquation& eqns, const dictionary& dict),
            (eqns, dict)
        );
        declareRunTimeSelectionTable
        (
            autoPtr,
            minimizationScheme,
            dictionaryMultivariate,
            (const scalarEquation& eqns, const dictionary& dict),
            (eqns, dict)
        );


    // Constructors

        //- Construct for given minimizationScheme
        minimizationScheme
        (
            const scalarEquation& eqns,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        minimizationScheme(const minimizationScheme&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<minimizationScheme> New
        (
            const scalarEquation& eqn,
            const dictionary& dict
        );


    //- Destructor
    virtual ~minimizationScheme()
    {}


    // Member Functions

        //- Return access to the tolerance
        inline scalarList& tolerances();

        //- Return the step number
        inline label nSteps() const;

        //- Return the errors
        inline scalarList& errors() const;

        //- Find the root of the equation
        virtual tmp<scalarField> solve() const;

        //- Find the root of the equation
        virtual tmp<scalarField> solve(const scalarField& x0) const;

        //- Find the root of the equation between xHigh and xLow
        virtual tmp<scalarField> solve
        (
            const scalarField& x0,
            const scalarField& xLow,
            const scalarField& xHigh
        ) const;

        //- Find the root of the equation
        virtual tmp<scalarField> solve
        (
            const scalarField& x0,
            const label li
        ) const;

        //- Find the root of the equation between xHigh and xLow
        virtual tmp<scalarField> solve
        (
            const scalarField& x0,
            const scalarField& xLow,
            const scalarField& xHigh,
            const label li
        ) const;

        //- Find the root of the equation between xHigh and xLow
        virtual tmp<scalarField> minimize
        (
            const scalarField& x0,
            const scalarField& xLow,
            const scalarField& xHigh,
            const label li
        ) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const minimizationScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "minimizationSchemeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
