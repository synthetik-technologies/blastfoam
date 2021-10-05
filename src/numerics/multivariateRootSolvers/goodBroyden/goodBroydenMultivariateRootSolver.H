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
    Foam::goodBroydenMultivariateRootSolver

Description
    Find multivariateMultivariateRoot using the goodBroyden method

SourceFiles
    multivariateMultivariateRootSolver.C

\*---------------------------------------------------------------------------*/

#ifndef goodBroydenMultivariateRootSolver_H
#define goodBroydenMultivariateRootSolver_H

#include "multivariateRootSolver.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class multivariateMultivariateRootSolver Declaration
\*---------------------------------------------------------------------------*/

class goodBroydenMultivariateRootSolver
:
    public multivariateRootSolver
{
    // Private data

        //- Finite difference value
        scalarField dx_;


public:

    //- Runtime type information
    TypeName("goodBroyden");

    // Constructors

        //- Construct for given multivariateMultivariateRootSolver
        goodBroydenMultivariateRootSolver
        (
            const scalarMultivariateEquation& eqn,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        goodBroydenMultivariateRootSolver
        (
            const multivariateRootSolver&
        ) = delete;


    //- Destructor
    virtual ~goodBroydenMultivariateRootSolver()
    {}


    // Member Functions

        //- Solve the multivariateMultivariateRoot
        virtual tmp<scalarField> findRoots
        (
            const scalarField& x0,
            const scalarField& xLow,
            const scalarField& xHigh,
            const label li
        ) const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const goodBroydenMultivariateRootSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
