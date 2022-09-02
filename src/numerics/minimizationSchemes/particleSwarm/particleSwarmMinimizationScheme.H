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
    Foam::particleSwarmMinimizationScheme

Description
    Gradient free minimization using a swarm of particles to randomly sample
    the parameter space.

References:
    \verbatim
        "Algorithms for Optimization",
        Kochenderfer, M. J., Wheeler, T. A.,
        The MIT Press, Cambridge, Massachusetts, 2019.
    \endverbatim

SourceFiles
    particleSwarmMinimizationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef particleSwarmMinimizationScheme_H
#define particleSwarmMinimizationScheme_H

#include "minimizationScheme.H"
#include "typeInfo.H"
#include "Random.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class particleSwarmMinimizationScheme Declaration
\*---------------------------------------------------------------------------*/

class particleSwarmMinimizationScheme
:
    public minimizationScheme
{
public:
    struct particle
    {
        scalarField x;
        scalar y;
        scalarField v;
        scalarField xBest;
        scalar yBest;
    };

private:

    //- Random number generator
    mutable Random rand_;

    //- Particles
    mutable List<particle> particles_;

    //- Coefficients
    scalar cLocal_;
    scalar cGlobal_;
    scalar vWeight_;


public:

    //- Runtime type information
    TypeName("particleSwarm");

    // Constructors

        //- Construct for given minimizationScheme
        particleSwarmMinimizationScheme
        (
            const scalarUnivariateEquation& eqn,
            const dictionary& dict
        );

    //- Destructor
    virtual ~particleSwarmMinimizationScheme()
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

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
