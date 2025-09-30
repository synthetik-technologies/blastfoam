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
    Foam::rootSolvers::univariate::Halley

Description
    Univariate root finding using the Halley method which requires first,
    and second derivatives.

SourceFiles
    HalleyUnivariateRootSolver.C

\*---------------------------------------------------------------------------*/

#ifndef HalleyUnivariateRootSolver_H
#define HalleyUnivariateRootSolver_H

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
                Class rootSolvers::univariate::Halley Declaration
\*---------------------------------------------------------------------------*/

class Halley
:
    public univariateRootSolver
{
public:

    //- Runtime type information
    TypeName("Halley");

    // Constructors

        //- Construct for given univariateRootSolver
        Halley
        (
            const scalarMultivariateEquation& eqn,
            const dictionary& dict
        );

        //- Copy constructor
        Halley
        (
            const scalarMultivariateEquation& eqn,
            const Halley& solver
        );

        //- Disallow default bitwise copy construction
        Halley(const Halley&) = delete;

        //- Clone the root solver
        virtual autoPtr<rootSolver> clone
        (
            const scalarMultivariateEquation& eqns
        ) const
        {
            return autoPtr<rootSolver>(new Halley(eqns, *this));
        }

        //- Clone the univariate root solver
        virtual autoPtr<univariateRootSolver> cloneUnivariate
        (
            const scalarMultivariateEquation& eqns
        ) const
        {
            return autoPtr<univariateRootSolver>(new Halley(eqns, *this));
        }


    //- Destructor
    virtual ~Halley();


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
