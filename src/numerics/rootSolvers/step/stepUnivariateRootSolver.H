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
    Foam::stepUnivariateRootSolver

Description
    Simple method to find a root of an equation by starting at the lower
    bound of an equation and stepping by the given distance until the the
    sign has changed. The step size is the reduced by the given factor.

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

    //- Reduction factor
    scalar f_;


public:

    //- Runtime type information
    TypeName("step");

    // Constructors

        //- Construct for given UnivariateRootSolver
        stepUnivariateRootSolver
        (
            const scalarMultivariateEquation& eqn,
            const dictionary& dict
        );

        //- Construct from equation and step size
        stepUnivariateRootSolver
        (
            const scalarMultivariateEquation& eqn,
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
