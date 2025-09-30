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
    Foam::rootSolvers::goodBroyden

Description
    Find multivariateRoot using the goodBroyden method. The inverse jacobian
    is explicitly update so that only one evaluation of the jacobian matrix
    and its inverse is needed

SourceFiles
    goodBroydenRootSolver.C

\*---------------------------------------------------------------------------*/

#ifndef goodBroydenRootSolver_H
#define goodBroydenRootSolver_H

#include "rootSolver.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace rootSolvers
{

/*---------------------------------------------------------------------------*\
                    Class rootSolvers::goodBroyden Declaration
\*---------------------------------------------------------------------------*/

class goodBroyden
:
    public rootSolver
{
    scalarField dX_;

public:

    //- Runtime type information
    TypeName("goodBroyden");

    // Constructors

        //- Construct for given multivariateRootSolver
        goodBroyden
        (
            const scalarMultivariateEquation& eqn,
            const dictionary& dict
        );

        //- Copy construct for given multivariateRootSolver
        goodBroyden
        (
            const scalarMultivariateEquation& eqn,
            const goodBroyden& solver
        );


        //- Disallow default bitwise copy construction
        goodBroyden
        (
            const goodBroyden&
        ) = delete;

        //- Clone the solver
        virtual autoPtr<rootSolver> clone
        (
            const scalarMultivariateEquation& eqns
        ) const
        {
            return autoPtr<rootSolver>(new goodBroyden(eqns, *this));
        }


    //- Destructor
    virtual ~goodBroyden();


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
