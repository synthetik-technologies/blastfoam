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
    Foam::minimizationSchemes::univariate::goldenRatio

Description
    Univariate minimization using the golden ratio to determine the next
    location

References:
    \verbatim
        "Algorithms for Optimization",
        Kochenderfer, M. J., Wheeler, T. A.,
        The MIT Press, Cambridge, Massachusetts, 2019.
    \endverbatim

SourceFiles
    goldenRatioUnivariateMinimizationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef goldenRatioUnivariateMinimizationScheme_H
#define goldenRatioUnivariateMinimizationScheme_H

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
        Class minimizationSchemes::univariate::goldenRatio Declaration
\*---------------------------------------------------------------------------*/

class goldenRatio
:
    public univariateMinimizationScheme
{
    const static scalar invPhi;
    const static scalar invPhi2;

public:

    //- Runtime type information
    TypeName("goldenRatio");

        //- Constructors
        goldenRatio
        (
            const scalarUnivariateEquation& eqn,
            const dictionary& dict
        );

        //- Copy constructor
        goldenRatio
        (
            const scalarUnivariateEquation& eqn,
            const goldenRatio& solver
        );

        //- Clone the root solver
        virtual autoPtr<minimizationScheme> clone
        (
            const scalarUnivariateEquation& eqns
        ) const
        {
            return autoPtr<minimizationScheme>(new goldenRatio(eqns, *this));
        }

        //- Clone the univariate root solver
        virtual autoPtr<univariateMinimizationScheme> cloneUnivariate
        (
            const scalarUnivariateEquation& eqns
        ) const
        {
            return autoPtr<univariateMinimizationScheme>
            (
                new goldenRatio(eqns, *this)
            );
        }


    //- Destructor
    virtual ~goldenRatio()
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
        void operator=(const goldenRatio&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace univariate
} // End namespace minimizationSchemes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
