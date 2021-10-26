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
    Foam::falsePointUnivariateRootSolver

Description
    Find root using the falsePoint method

SourceFiles
    falsePointUnivariateRootSolver.C

\*---------------------------------------------------------------------------*/

#ifndef falsePointUnivariateRootSolver_H
#define falsePointUnivariateRootSolver_H

#include "univariateRootSolver.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class falsePointUnivariateRootSolver Declaration
\*---------------------------------------------------------------------------*/

class falsePointUnivariateRootSolver
:
    public univariateRootSolver
{
public:

    //- Runtime type information
    TypeName("falsePoint");

    // Constructors

        //- Construct for given univariateRootSolver
        falsePointUnivariateRootSolver
        (
            const multivariateEquation<scalar>& eqn,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        falsePointUnivariateRootSolver(const univariateRootSolver&) = delete;


    //- Destructor
    virtual ~falsePointUnivariateRootSolver();


    // Member Functions

        //- Solve the root
        virtual scalar findRoot
        (
            const scalar x0,
            const scalar xLow,
            const scalar xHigh,
            const label li
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
