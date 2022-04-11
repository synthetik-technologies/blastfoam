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
    Foam::univariateMinimizationScheme

Description
    General class for minimizing univariate, scalar equations

SourceFiles
    univariateMinimizationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef univariateMinimizationScheme_H
#define univariateMinimizationScheme_H

#include "EquationsFwd.H"
#include "minimizationScheme.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class univariateMinimizationScheme Declaration
\*---------------------------------------------------------------------------*/

class univariateMinimizationScheme
:
    public minimizationScheme
{

protected:

    // Protected data

        //- Reference to scalarEquation
        const scalarEquation& eqn_;

        //- Number of location to pre sample
        label nSample_;

        //- Is the solution converged in x
        bool convergedXScale(const scalar error, const scalar s) const;

        //- Is the solution converged in x
        bool convergedX(const scalar x1, const scalar x2) const;

        //- Is the solution converged in y
        bool convergedYScale(const scalar error, const scalar s) const;

        //- Is the solution converged in y
        bool convergedY(const scalar y1, const scalar y2) const;

        //- Print information about the current step
        void printStepInformation(const scalar val) const;

        //- Print information about the current step
        scalar printFinalInformation(const scalar val) const;

        //- Sample and return the interval with the smallest value
        void sample(scalar&, scalar&, const label) const;

public:

    //- Runtime type information
    TypeNameNoDebug("univariateMinimizationScheme");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            univariateMinimizationScheme,
            dictionaryZero,
            (const scalarUnivariateEquation& eqn, const dictionary& dict),
            (eqn, dict)
        );
        declareRunTimeSelectionTable
        (
            autoPtr,
            univariateMinimizationScheme,
            dictionaryOne,
            (const scalarUnivariateEquation& eqn, const dictionary& dict),
            (eqn, dict)
        );
        declareRunTimeSelectionTable
        (
            autoPtr,
            univariateMinimizationScheme,
            dictionaryTwo,
            (const scalarUnivariateEquation& eqn, const dictionary& dict),
            (eqn, dict)
        );


    // Constructors

        //- Construct for given univariateMinimizationScheme
        univariateMinimizationScheme
        (
            const scalarUnivariateEquation& eqn,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        univariateMinimizationScheme(const univariateMinimizationScheme&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<univariateMinimizationScheme> New
        (
            const scalarUnivariateEquation& eqn,
            const dictionary& dict
        );


    //- Destructor
    virtual ~univariateMinimizationScheme()
    {}


    // Member Functions

        //- Return the x tolerance
        inline scalar xTolerance() const;

        //- Return the y tolerance
        inline scalar yTolerance() const;

        //- Return the x error
        inline scalar xError() const;

        //- Return the y error
        inline scalar yError() const;

        using minimizationScheme::solve;

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

        //- Find the root of the equation between xHigh and xLow
        virtual tmp<scalarField> minimize
        (
            const scalarList& x0,
            const scalarList& xLow,
            const scalarList& xHigh,
            const label li
        ) const;


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
