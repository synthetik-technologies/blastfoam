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
namespace rootSolvers
{
namespace univariate
{
/*---------------------------------------------------------------------------*\
             Class rootSolvers::univariate::NewtonRaphson Declaration
\*---------------------------------------------------------------------------*/

class NewtonRaphson
:
    public univariateRootSolver
{
public:

    //- Runtime type information
    TypeName("NewtonRaphson");

    // Constructors

        //- Construct for given univariateRootSolver
        NewtonRaphson
        (
            const scalarMultivariateEquation& eqn,
            const dictionary& dict
        );

        //- Copy constructor
        NewtonRaphson
        (
            const scalarMultivariateEquation& eqn,
            const NewtonRaphson& solver
        );

        //- Disallow default bitwise copy construction
        NewtonRaphson(const NewtonRaphson&) = delete;

        //- Clone the root solver
        virtual autoPtr<rootSolver> clone
        (
            const scalarMultivariateEquation& eqns
        ) const
        {
            return autoPtr<rootSolver>(new NewtonRaphson(eqns, *this));
        }

        //- Clone the univariate root solver
        virtual autoPtr<univariateRootSolver> cloneUnivariate
        (
            const scalarMultivariateEquation& eqns
        ) const
        {
            return autoPtr<univariateRootSolver>
            (
                new NewtonRaphson(eqns, *this)
            );
        }


    //- Destructor
    virtual ~NewtonRaphson();


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

} // End namespace univariate
} // End namespace rootSolvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
