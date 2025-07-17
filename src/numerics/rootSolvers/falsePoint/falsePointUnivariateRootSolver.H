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
    Foam::rootSolvers::univariate::falsePoint

Description
    Univariate root finding using the false-point method.

SourceFiles
    falsePointUnivariateRootSolver.C

\*---------------------------------------------------------------------------*/

#ifndef falsePointUnivariateRootSolver_H
#define falsePointUnivariateRootSolver_H

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
                Class rootSolvers::univariate::falsePoint Declaration
\*---------------------------------------------------------------------------*/

class falsePoint
:
    public univariateRootSolver
{
public:

    //- Runtime type information
    TypeName("falsePoint");

    // Constructors

        //- Construct for given univariateRootSolver
        falsePoint
        (
            const scalarMultivariateEquation& eqn,
            const dictionary& dict
        );

        //- Copy constructor
        falsePoint
        (
            const scalarMultivariateEquation& eqn,
            const falsePoint& solver
        );

        //- Disallow default bitwise copy construction
        falsePoint(const falsePoint&) = delete;

        //- Clone the root solver
        virtual autoPtr<rootSolver> clone
        (
            const scalarMultivariateEquation& eqns
        ) const
        {
            return autoPtr<rootSolver>(new falsePoint(eqns, *this));
        }

        //- Clone the univariate root solver
        virtual autoPtr<univariateRootSolver> cloneUnivariate
        (
            const scalarMultivariateEquation& eqns
        ) const
        {
            return autoPtr<univariateRootSolver>(new falsePoint(eqns, *this));
        }


    //- Destructor
    virtual ~falsePoint();


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
} // End name space rootSolvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
