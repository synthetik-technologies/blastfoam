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
    Foam::minimizationSchemes::univariate::quadraticFit

Description
    Univariate minimization using an imposed quadratic fit

References:
    \verbatim
        "Algorithms for Optimization",
        Kochenderfer, M. J., Wheeler, T. A.,
        The MIT Press, Cambridge, Massachusetts, 2019.
    \endverbatim

SourceFiles
    quadraticFitUnivariateMinimizationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef quadraticFitUnivariateMinimizationScheme_H
#define quadraticFitUnivariateMinimizationScheme_H

#include "univariateMinimizationScheme.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace minimizationSchemes
{
namespace univariate
{

/*---------------------------------------------------------------------------*\
                        Class quadraticFitUnivariateMinimizationScheme Declaration
\*---------------------------------------------------------------------------*/

class quadraticFit
:
    public univariateMinimizationScheme
{
public:

    //- Runtime type information
    TypeName("quadraticFit");

    // Constructors

        //- Construct for given univariateMinimizationScheme
        quadraticFit
        (
            const scalarUnivariateEquation& eqn,
            const dictionary& dict
        );

        //- Copy constructor
        quadraticFit
        (
            const scalarUnivariateEquation& eqn,
            const quadraticFit& solver
        );

        //- Clone the root solver
        virtual autoPtr<minimizationScheme> clone
        (
            const scalarUnivariateEquation& eqns
        ) const
        {
            return autoPtr<minimizationScheme>(new quadraticFit(eqns, *this));
        }

        //- Clone the univariate root solver
        virtual autoPtr<univariateMinimizationScheme> cloneUnivariate
        (
            const scalarUnivariateEquation& eqns
        ) const
        {
            return autoPtr<univariateMinimizationScheme>
            (
                new quadraticFit(eqns, *this)
            );
        }


    //- Destructor
    virtual ~quadraticFit()
    {}


    // Member Functions

        //- Solve the root
        virtual scalar minimize
        (
            const scalar x0,
            const scalar x1,
            const scalar x2,
            const label li
        ) const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const quadraticFit&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace univariate
} // End namespace minimizationSchemes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
