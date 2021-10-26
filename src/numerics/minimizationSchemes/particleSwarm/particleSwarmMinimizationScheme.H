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
    Foam::particleSwarmMinimizationScheme

Description
    Find root using the particleSwarm method

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

    //- Name of the error
    virtual word errorName() const
    {
        return "variance";
    }

public:

    //- Runtime type information
    TypeName("particleSwarm");

    // Constructors

        //- Construct for given minimizationScheme
        particleSwarmMinimizationScheme
        (
            const scalarEquation& eqn,
            const dictionary& dict
        );

    //- Destructor
    virtual ~particleSwarmMinimizationScheme()
    {}


    // Member Functions

        //- Solve the root
        virtual tmp<scalarField> minimize
        (
            const scalarField& x0,
            const scalarField& x1,
            const scalarField& x2,
            const label li
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
