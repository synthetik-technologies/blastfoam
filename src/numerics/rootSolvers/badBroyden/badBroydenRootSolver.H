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
    Foam::rootSolvers::badBroyden

Description
    Find multivariateRoot using the badBroyden method. The inverse jacobian
    is explicitly update so that only one evaluation of the jacobian matrix
    and its inverse is needed

SourceFiles
    badBroydenRootSolver.C

\*---------------------------------------------------------------------------*/

#ifndef badBroydenRootSolver_H
#define badBroydenRootSolver_H

#include "rootSolver.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace rootSolvers
{

/*---------------------------------------------------------------------------*\
                          Class rootSolvers::badBroyden Declaration
\*---------------------------------------------------------------------------*/

class badBroyden
:
    public rootSolver
{
    scalarField dX_;


public:

    //- Runtime type information
    TypeName("badBroyden");

    // Constructors

        //- Construct for given multivariateRootSolver
        badBroyden
        (
            const scalarMultivariateEquation& eqn,
            const dictionary& dict
        );

        //- Copy construct for given multivariateRootSolver
        badBroyden
        (
            const scalarMultivariateEquation& eqn,
            const badBroyden& solver
        );

        //- Disallow default bitwise copy construction
        badBroyden(const badBroyden&) = delete;

        //- Clone the solver
        virtual autoPtr<rootSolver> clone
        (
            const scalarMultivariateEquation& eqns
        ) const
        {
            return autoPtr<rootSolver>(new badBroyden(eqns, *this));
        }



    //- Destructor
    virtual ~badBroyden();


    // Member Functions

        //- Solve the multivariateRoot
        virtual tmp<scalarField> findRoots
        (
            const scalarList& x0,
            const scalarList& xLow,
            const scalarList& xHigh,
            const label li
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace rootSolvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
