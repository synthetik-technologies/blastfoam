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
    Foam::NewtonRaphsonRootSolver

Description
    Multivariate root finding using a Newton-Raphson method using the
    Jacobian of the system

SourceFiles
    NewtonRaphsonRootSolver.C

\*---------------------------------------------------------------------------*/

#ifndef NewtonRaphsonRootSolver_H
#define NewtonRaphsonRootSolver_H

#include "rootSolver.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class NewtonRaphsonRootSolver Declaration
\*---------------------------------------------------------------------------*/

class NewtonRaphsonRootSolver
:
    public rootSolver
{
// Private data

    //- Damping factor
    scalar beta_;


public:

    //- Runtime type information
    TypeName("NewtonRaphson");

    // Constructors

        //- Construct for given multivariateMultivariateRootSolver
        NewtonRaphsonRootSolver
        (
            const scalarMultivariateEquation& eqn,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        NewtonRaphsonRootSolver
        (
            const NewtonRaphsonRootSolver&
        ) = delete;


    //- Destructor
    virtual ~NewtonRaphsonRootSolver();


    // Member Functions

        //- Solve the multivariateMultivariateRoot
        virtual tmp<scalarField> findRoots
        (
            const scalarField& x0,
            const scalarField& xLow,
            const scalarField& xHigh,
            const label li
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
