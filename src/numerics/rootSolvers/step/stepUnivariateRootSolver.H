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
    Foam::rootSolvers::univariate::step

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
namespace rootSolvers
{
namespace univariate
{

/*---------------------------------------------------------------------------*\
                Class rootSolvers::univariate::step Declaration
\*---------------------------------------------------------------------------*/

class step
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
        step
        (
            const scalarMultivariateEquation& eqn,
            const dictionary& dict
        );

        //- Construct from equation and step size
        step
        (
            const scalarMultivariateEquation& eqn,
            const scalar dx
        );

        //- Copy constructor
        step
        (
            const scalarMultivariateEquation& eqn,
            const step& solver
        );

        //- Disallow default bitwise copy construction
        step(const step&) = delete;

        //- Clone the root solver
        virtual autoPtr<rootSolver> clone
        (
            const scalarMultivariateEquation& eqns
        ) const
        {
            return autoPtr<rootSolver>(new step(eqns, *this));
        }

        //- Clone the univariate root solver
        virtual autoPtr<univariateRootSolver> cloneUnivariate
        (
            const scalarMultivariateEquation& eqns
        ) const
        {
            return autoPtr<univariateRootSolver>(new step(eqns, *this));
        }


    //- Destructor
    virtual ~step();


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
