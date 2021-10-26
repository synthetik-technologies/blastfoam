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
    Foam::badBroydenRootSolver

Description
    Find multivariateRoot using the badBroyden method. The inverse jacobian
    is explicitly update so that only one evaluation of the jacobian matrix
    and its inverse is needed

SourceFiles
    multivariateRootSolver.C

\*---------------------------------------------------------------------------*/

#ifndef badBroydenRootSolver_H
#define badBroydenRootSolver_H

#include "rootSolver.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class multivariateRootSolver Declaration
\*---------------------------------------------------------------------------*/

class badBroydenRootSolver
:
    public rootSolver
{
public:

    //- Runtime type information
    TypeName("badBroyden");

    // Constructors

        //- Construct for given multivariateRootSolver
        badBroydenRootSolver
        (
            const multivariateEquation<scalar>& eqn,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        badBroydenRootSolver(const badBroydenRootSolver&) = delete;


    //- Destructor
    virtual ~badBroydenRootSolver();


    // Member Functions

        //- Solve the multivariateRoot
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
