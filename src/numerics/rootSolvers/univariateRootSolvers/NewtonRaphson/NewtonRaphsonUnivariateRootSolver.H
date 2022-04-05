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
    Foam::NewtonRaphsonUnivariateRootSolver

Description
    Univariate root finding using a NewtonRaphson method. Requires the first
    derivative of an equation.

SourceFiles
    NewtonRaphsonUnivariateRootSolver.C

\*---------------------------------------------------------------------------*/

#ifndef NewtonRaphsonUnivariateRootSolver_H
#define NewtonRaphsonUnivariateRootSolver_H

#include "univariateRootSolver.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class NewtonRaphsonUnivariateUnivariateRootSolver Declaration
\*---------------------------------------------------------------------------*/

class NewtonRaphsonUnivariateRootSolver
:
    public univariateRootSolver
{
public:

    //- Runtime type information
    TypeName("NewtonRaphson");

    // Constructors

        //- Construct for given univariateRootSolver
        NewtonRaphsonUnivariateRootSolver
        (
            const scalarMultivariateEquation& eqn,
            const dictionary& dict
        );

        //- Construct for given univariateRootSolver
        NewtonRaphsonUnivariateRootSolver
        (
            const scalarMultivariateEquation& eqn,
            const scalar tolerance = 1e-6,
            const label maxSteps = 100
        );

        //- Disallow default bitwise copy construction
        NewtonRaphsonUnivariateRootSolver(const univariateRootSolver&) = delete;


    //- Destructor
    virtual ~NewtonRaphsonUnivariateRootSolver();


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
