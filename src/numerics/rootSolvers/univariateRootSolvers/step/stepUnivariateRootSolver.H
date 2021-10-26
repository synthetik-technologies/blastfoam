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
    Foam::stepUnivariateRootSolver

Description
    Find root using the step method

SourceFiles
    stepUnivariateRootSolver.C

\*---------------------------------------------------------------------------*/

#ifndef stepUnivariateRootSolver_H
#define stepUnivariateRootSolver_H

#include "univariateRootSolver.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class stepRootSolver Declaration
\*---------------------------------------------------------------------------*/

class stepUnivariateRootSolver
:
    public univariateRootSolver
{
// Private data

    //- Initial step size
    scalar dx_;


public:

    //- Runtime type information
    TypeName("step");

    // Constructors

        //- Construct for given UnivariateRootSolver
        stepUnivariateRootSolver
        (
            const multivariateEquation<scalar>& eqn,
            const dictionary& dict
        );

        //- Construct from equation and step size
        stepUnivariateRootSolver
        (
            const equation& eqn,
            const scalar dx
        );


        //- Disallow default bitwise copy construction
        stepUnivariateRootSolver(const rootSolver&) = delete;


    //- Destructor
    virtual ~stepUnivariateRootSolver();


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
