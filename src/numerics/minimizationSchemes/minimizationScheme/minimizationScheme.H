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
    Foam::minimizationScheme

Description
    General class for minimizing multivariate, scalar equations

SourceFiles
    minimizationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef minimizationScheme_H
#define minimizationScheme_H

#include "UnivariateEquationsFwd.H"
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
        const scalarUnivariateEquation& eqns_;

        //- Absolute convergence tolerances for x
        scalarField xTolerances_;

        //- Relative convergence tolerances for x
        scalarField xRelTolerances_;

        //- Absolute convergence tolerances for y
        scalarField yTolerances_;

        //- Relative convergence tolerances for y
        scalarField yRelTolerances_;

        //- The maximum number of sub-steps allowed for the integration step
        label maxSteps_;

        //- Current step number
        mutable label stepi_;

        //- Current absolute x errors
        mutable scalarList xErrors_;

        //- Current relative x errors
        mutable scalarList xRelErrors_;

        //- Current absolute delta y errors
        mutable scalarList yErrors_;

        //- Current relative delta y errors
        mutable scalarList yRelErrors_;

        //- Number of pre-samples in each direction
        labelList nSamples_;

        //- Is the error normalized
        Switch normalize_;

        //- Line search reduction coefficient
        scalar tau_;

        //- In convergence in the delta checked
        bool checkY_;


        //- Initial line search parameter
        virtual scalar alpha
        (
            const scalarList& grad,
            const scalarList& gradOld
        ) const;

        //- Line search
        scalar lineSearch
        (
            const scalarList& x0,
            const scalarList& grad,
            const scalarList& gradOld,
            const label li,
            scalar& fx
        ) const;

        //- Line search
        scalar lineSearch
        (
            const scalarList& x0,
            const scalarList& grad,
            const label li,
            scalar& fx
        ) const;

        //- Return the L2 norm the given list
        scalar norm(const scalarList& lst) const;

        //- Return lst^T * lst
        scalar inner(const scalarList& lst1, const scalarList& lst2) const;

        //- Is the solution converged in w.r.t x
        bool convergedXScale
        (
            const scalarList& error,
            const scalarList& s
        ) const;

        bool converged
        (
            const scalarList& absErrors,
            const scalarList& relErrors,
            const scalarList& absTols,
            const scalarList& relTols
        ) const;

        //- Is the solution converged in w.r.t x
        bool convergedX(const scalarList& x1, const scalarList& x2) const;

        //- Is the solution converged in w.r.t y
        bool convergedYScale
        (
            const scalarList& error,
            const scalarList& s
        ) const;

        //- Is the solution converged in w.r.t y
        bool convergedY(const scalarList& y1, const scalarList& y2) const;

        //- Print information about the current step
        void printStepInformation(const scalarList& vals) const;

        //- Print information about the current step
        void printFinalInformation(const scalarList& vals) const;

        //- Sample and return the interval with the smallest value
        void sample
        (
            scalarList&,
            scalarList&,
            labelList& xBest,
            scalar& yBest,
            labelList& indicies,
            const label,
            const label diri
        ) const;


public:

    //- Runtime type information
    TypeName("minimizationScheme");

    // Declare run-time constructor selection tables
        declareRunTimeSelectionTable
        (
            autoPtr,
            minimizationScheme,
            dictionaryUnivariate,
            (const scalarUnivariateEquation& eqns, const dictionary& dict),
            (eqns, dict)
        );
        declareRunTimeSelectionTable
        (
            autoPtr,
            minimizationScheme,
            dictionaryMultivariate,
            (const scalarUnivariateEquation& eqns, const dictionary& dict),
            (eqns, dict)
        );


    // Constructors

        //- Construct for given minimizationScheme
        minimizationScheme
        (
            const scalarUnivariateEquation& eqns,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        minimizationScheme(const minimizationScheme&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<minimizationScheme> New
        (
            const scalarUnivariateEquation& eqn,
            const dictionary& dict
        );


    //- Destructor
    virtual ~minimizationScheme()
    {}


    // Member Functions

        //- Return the step number
        inline label nSteps() const;

        //- Return the absolute x errors
        inline scalarList& xErrors() const;

        //- Return the relative x errors
        inline scalarList& xRelErrors() const;

        //- Return the absolute y errors
        inline scalarList& yErrors() const;

        //- Return the relative y errors
        inline scalarList& yRelErrors() const;

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
        void operator=(const minimizationScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "minimizationSchemeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
