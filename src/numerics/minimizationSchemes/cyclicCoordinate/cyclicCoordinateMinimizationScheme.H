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
    Foam::minimizationSchemes::cyclicCoordinate

Description


References:
    \verbatim
        "Algorithms for Optimization",
        Kochenderfer, M. J., Wheeler, T. A.,
        The MIT Press, Cambridge, Massachusetts, 2019.
    \endverbatim

SourceFiles
    cyclicCoordinateMinimizationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicCoordinateMinimizationScheme_H
#define cyclicCoordinateMinimizationScheme_H

#include "basisMinimizationScheme.H"
#include "cmptLineSearchEquation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace minimizationSchemes
{

/*---------------------------------------------------------------------------*\
                    Class cyclicCoordinate Declaration
\*---------------------------------------------------------------------------*/

class cyclicCoordinate
:
    public basis
{
    // Private member data

        //- Use acceleration
        bool accelerate_;


public:

    //- Runtime type information
    TypeName("cyclicCoordinate");

    // Constructors

        //- Construct for given minimizationScheme
        cyclicCoordinate
        (
            const scalarUnivariateEquation& eqn,
            const dictionary& dict
        );

        //- Copy constructor
        cyclicCoordinate
        (
            const scalarUnivariateEquation& eqn,
            const cyclicCoordinate& solver
        );

        //- Clone the root solver
        virtual autoPtr<minimizationScheme> clone
        (
            const scalarUnivariateEquation& eqns
        ) const
        {
            return autoPtr<minimizationScheme>
            (
                new cyclicCoordinate(eqns, *this)
            );
        }


    //- Destructor
    virtual ~cyclicCoordinate()
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
