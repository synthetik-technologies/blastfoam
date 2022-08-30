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
    Foam::bisectionUnivariateMinimizationScheme

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

/*---------------------------------------------------------------------------*\
                        Class bisectionUnivariateMinimizationScheme Declaration
\*---------------------------------------------------------------------------*/

class bisectionUnivariateMinimizationScheme
:
    public univariateMinimizationScheme
{
public:

    //- Runtime type information
    TypeName("bisection");

    // Constructors

        //- Construct for given univariateUnivariateMinimizationScheme
        bisectionUnivariateMinimizationScheme
        (
            const scalarUnivariateEquation& eqn,
            const dictionary& dict
        );


    //- Destructor
    virtual ~bisectionUnivariateMinimizationScheme()
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
        void operator=(const bisectionUnivariateMinimizationScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
