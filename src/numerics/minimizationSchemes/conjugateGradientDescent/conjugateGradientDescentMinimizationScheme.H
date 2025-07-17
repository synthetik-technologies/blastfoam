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
    Foam::minimizationSchemes::conjugateGradientDescent

Description
    Minimization scheme using either gradient computed using a finite
    gradient of an implemented, analytical version

References:
    \verbatim
        "Algorithms for Optimization",
        Kochenderfer, M. J., Wheeler, T. A.,
        The MIT Press, Cambridge, Massachusetts, 2019.
    \endverbatim

SourceFiles
    conjugateGradientDescentMinimizationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef conjugateGradientDescentMinimizationScheme_H
#define conjugateGradientDescentMinimizationScheme_H

#include "minimizationScheme.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace minimizationSchemes
{

/*---------------------------------------------------------------------------*\
                Class conjugateGradientDescent Declaration
\*---------------------------------------------------------------------------*/

class conjugateGradientDescent
:
    public minimizationScheme
{
public:

    //- Runtime type information
    TypeName("conjugateGradientDescent");

    // Constructors

        //- Construct for given minimizationScheme
        conjugateGradientDescent
        (
            const scalarUnivariateEquation& eqn,
            const dictionary& dict
        );

        //- Copy constructor
        conjugateGradientDescent
        (
            const scalarUnivariateEquation& eqn,
            const conjugateGradientDescent& solver
        );

        //- Clone the root solver
        virtual autoPtr<minimizationScheme> clone
        (
            const scalarMultivariateEquation& eqns
        ) const
        {
            return autoPtr<minimizationScheme>
            (
                new conjugateGradientDescent(eqns, *this)
            );
        }


    //- Destructor
    virtual ~conjugateGradientDescent()
    {}


    // Member Functions

        //- Solve the root
        virtual tmp<scalarField> minimize
        (
            const scalarList& x0,
            const scalarList& x1,
            const scalarList& x2,
            const label li
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace minimizationSchemes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
