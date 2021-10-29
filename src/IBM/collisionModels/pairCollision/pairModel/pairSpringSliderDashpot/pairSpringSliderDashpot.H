/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::pairSpringSliderDashpot

Description
    Pair forces between particles colliding with a spring, slider, damper model

\*---------------------------------------------------------------------------*/

#ifndef immersedPairSpringSliderDashpot_H
#define immersedPairSpringSliderDashpot_H

#include "pairModel.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class pairSpringSliderDashpot Declaration
\*---------------------------------------------------------------------------*/

class pairSpringSliderDashpot
:
    public pairModel
{
    // Private Data

        //- Effective Young's modulus value, assuming both particles have
        //  the same E value
        scalar Estar_;

        //- Effective shear modulus value, assuming both particles have
        //  the same Poisson's ratio and Young's modulus
        scalar Gstar_;

        //- alpha-coefficient, related to coefficient of restitution
        scalar alpha_;

        //- Spring power (b = 1 for linear, b = 3/2 for Hertzian)
        scalar b_;

        //- Coefficient of friction in for tangential sliding
        scalar mu_;

        //- Cohesion energy density [J/m^3]
        scalar cohesionEnergyDensity_;

        // Switch cohesion on and off
        bool cohesion_;

        //- The number of steps over which to resolve the minimum
        //  harmonic approximation of the collision period
        scalar collisionResolutionSteps_;

        //- Volume factor for determining the equivalent size of a
        //  parcel where nParticles is not 1.  The equivalent size of
        //  the parcel is
        //      parcelEquivVolume = volumeFactor*nParticles*p.volume()
        //  so
        //      parcelEquivD = cbrt(volumeFactor*nParticles)*p.d()
        //  + When volumeFactor = 1, the particles are compressed
        //    together so that the equivalent volume of the parcel is
        //    the sum of the constitutent particles
        //  + When volumeFactor = 3*sqrt(2)/pi, the particles are
        //    close packed, but uncompressed.
        //  + When volumeFactor > 3*sqrt(2)/pi, the particles loosely
        //    grouped.
        // 3*sqrt(2)/pi = 1.350474 is the volume factor for close
        // packing, i.e pi/(3*sqrt(2)) is the maximum close packing
        // factor
        scalar volumeFactor_;

        //- Switch to control use of equivalent size particles.  Used
        //  because the calculation can be very expensive.
        bool useEquivalentSize_;


public:

    //- Runtime type information
    TypeName("pairSpringSliderDashpot");


    // Constructors

        //- Construct from dictionary
        pairSpringSliderDashpot(const dictionary& dict);


    //- Destructor
    virtual ~pairSpringSliderDashpot();


    // Member Functions

        //- Return the volumeFactor
        inline scalar volumeFactor() const
        {
            return volumeFactor_;
        }

        // Return the area of overlap between two spheres of radii rA and rB,
        // centres separated by a distance rAB.  Assumes rAB < (rA + rB).
        inline scalar overlapArea(scalar rA, scalar rB, scalar rAB) const
        {
            // Reference:
            // http://mathworld.wolfram.com/Sphere-SphereIntersection.html
            return
                constant::mathematical::pi/4.0
               /sqr(rAB)
               *(
                    (-rAB + rA - rB)
                   *(-rAB - rA + rB)
                   *(-rAB + rA + rB)
                   *( rAB + rA + rB)
                );
        }

        //- Calculate the pair interaction between parcels
        virtual void evaluatePair
        (
            immersedBoundaryObject& OA,
            immersedBoundaryObject& OB,
            const pairCollisionData& map
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
