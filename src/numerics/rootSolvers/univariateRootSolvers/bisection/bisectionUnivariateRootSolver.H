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
    Foam::bisectionUnivariateRootSolver

Description
    Find root using the bisection method

SourceFiles
    bisectionUnivariateRootSolver.C

\*---------------------------------------------------------------------------*/

#ifndef bisectionUnivariateRootSolver_H
#define bisectionUnivariateRootSolver_H

#include "univariateRootSolver.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class bisectionUnivariateRootSolver Declaration
\*---------------------------------------------------------------------------*/

class bisectionUnivariateRootSolver
:
    public univariateRootSolver
{
public:

    //- Runtime type information
    TypeName("bisection");

    // Constructors

        //- Construct for given univariateRootSolver
        bisectionUnivariateRootSolver
        (
            const multivariateEquation<scalar>& eqn,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        bisectionUnivariateRootSolver(const univariateRootSolver&) = delete;


    //- Destructor
    virtual ~bisectionUnivariateRootSolver();


    // Member Functions

        //- Solve the root
        virtual scalar findRoot
        (
            const scalar x0,
            const scalar x1,
            const scalar x2,
            const label li
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
