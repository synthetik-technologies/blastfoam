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
    Foam::wallSpringSliderDashpot

Description
    Forces between particles and walls, interacting with a spring,
    slider, damper model

\*---------------------------------------------------------------------------*/

#ifndef immersedWallSpringSliderDashpot_H
#define immersedWallSpringSliderDashpot_H

#include "wallModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class wallSpringSliderDashpot Declaration
\*---------------------------------------------------------------------------*/

class wallSpringSliderDashpot
:
    public wallModel
{
    // Private Data

        //- Effective Young's modulus value
        scalar Estar_;

        //- Effective shear modulus value
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


    // Private Member Functions

        //- Calculate the wall interaction for a parcel at a given site
//         void evaluateWall
//         (
//             typename CloudType::parcelType& p,
//             const point& site,
//             const WallSiteData<vector>& data,
//             scalar pREff,
//             scalar kN,
//             bool cohesion
//         ) const;


public:

    //- Runtime type information
    TypeName("wallSpringSliderDashpot");


    // Constructors

        //- Construct from dictionary
        wallSpringSliderDashpot(const dictionary& dict);


    //- Destructor
    virtual ~wallSpringSliderDashpot();


    // Member Functions

        //- Return the volumeFactor
        inline scalar volumeFactor() const
        {
            return volumeFactor_;
        }

        //- Calculate the wall interaction for a parcel
        virtual void evaluateWall
        (
            const polyPatch& patch,
            immersedBoundaryObject& object,
            const wallCollisionData& map
        ) const;

        //- Calculate the wall interaction for a parcel
        virtual void evaluateExternalForce
        (
            const polyPatch& patch,
            immersedBoundaryObject& object,
            const wallCollisionData& map
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
