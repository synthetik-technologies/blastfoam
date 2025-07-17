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
    Foam::univariateMinimizationSchemes::bisection

Description
    Univariate minimization using the bisection method

SourceFiles
    bisectionUnivariateMinimizationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef bisectionUnivariateMinimizationScheme_H
#define bisectionUnivariateMinimizationScheme_H

#include "univariateMinimizationScheme.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace minimizationSchemes
{
namespace univariate
{

/*---------------------------------------------------------------------------*\
                    Class bisection Declaration
\*---------------------------------------------------------------------------*/

class bisection
:
    public univariateMinimizationScheme
{
public:

    //- Runtime type information
    TypeName("bisection");

    // Constructors

        //- Constructor
        bisection
        (
            const scalarUnivariateEquation& eqn,
            const dictionary& dict
        );

        //- Copy constructor
        bisection
        (
            const scalarUnivariateEquation& eqn,
            const bisection& solver
        );

        //- Clone the root solver
        virtual autoPtr<minimizationScheme> clone
        (
            const scalarUnivariateEquation& eqns
        ) const
        {
            return autoPtr<minimizationScheme>(new bisection(eqns, *this));
        }

        //- Clone the univariate root solver
        virtual autoPtr<univariateMinimizationScheme> cloneUnivariate
        (
            const scalarUnivariateEquation& eqns
        ) const
        {
            return autoPtr<univariateMinimizationScheme>
            (
                new bisection(eqns, *this)
            );
        }


    //- Destructor
    virtual ~bisection()
    {}


    // Member Functions

        //- Solve the root
        virtual scalar minimize
        (
            const scalar x0,
            const scalar x1,
            const scalar x2,
            const label li
        ) const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const bisection&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace univariate
} // End namespace minimizationSchemes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
