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

\*---------------------------------------------------------------------------*/

#include "pairSpringSliderDashpot.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pairSpringSliderDashpot, 0);
    addToRunTimeSelectionTable
    (
        pairModel,
        pairSpringSliderDashpot,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pairSpringSliderDashpot::pairSpringSliderDashpot
(
    const dictionary& dict
)
:
    pairModel(dict, typeName),
    Estar_(),
    Gstar_(),
    alpha_(this->coeffDict().template lookup<scalar>("alpha")),
    b_(this->coeffDict().template lookup<scalar>("b")),
    mu_(this->coeffDict().template lookup<scalar>("mu")),
    cohesionEnergyDensity_
    (
        this->coeffDict().template lookup<scalar>("cohesionEnergyDensity")
    ),
    cohesion_(false)
{

    scalar nu = this->coeffDict().template lookup<scalar>("nu");

    scalar E = this->coeffDict().template lookup<scalar>("E");

    Estar_ = E/(2.0*(1.0 - sqr(nu)));

    scalar G = E/(2.0*(1.0 + nu));

    Gstar_ = G/(2.0*(2.0 - nu));

    cohesion_ = (mag(cohesionEnergyDensity_) > vSmall);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pairSpringSliderDashpot::~pairSpringSliderDashpot()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::pairSpringSliderDashpot::evaluatePair
(
    immersedBoundaryObject& objectA,
    immersedBoundaryObject& objectB,
    const pairCollisionData& map
) const
{
    scalar dT = objectA.pMesh().time().deltaTValue();
    scalar massA = objectA.mass();
    scalar massB = objectB.mass();

    tensor invIA = objectA.invMomentOfInertia();
    tensor invIB = objectB.invMomentOfInertia();

    vector xcA = objectA.centreOfRotation();
    vector xcB = objectB.centreOfRotation();

    forAll(map, pairi)
    {
        const label facei = map[pairi][0];
        const label facej = map[pairi][1];
        if (facei < 0 || facej < 0)
        {
            continue;
        }

        const vector& xA = objectA.faceCentres()[facei];
        const vector& xB = objectB.faceCentres()[facej];

        vector hitPoint((xA + xB)*0.5);
        vector normal(xcB - xcA);
        normal /= mag(normal) + vSmall;

        scalar delta = (xB - xA) & normal;

        if (delta < 0)
        {
            delta = mag(delta);

            // Particles in collision
            vector vAB(objectA.velocity(hitPoint) - objectB.velocity(hitPoint));

            // Effective mass
            vector rA(xA - xcA);
            vector rB(xB - xcB);
            scalar M =
                1.0
               /(
                    1.0/massA + 1.0/massB
                  + (normal & (invIA & ((rA ^ normal) ^ rA)))
                  + (normal & (invIB & ((rB ^ normal) ^ rB)))
                );

            scalar areaA(mag(objectA.Sf()[facei]));
            scalar areaB(mag(objectB.Sf()[facej]));
            scalar LAB
            (
                0.5*(areaA/(mag(rA) + vSmall) + areaB/(mag(rB) + vSmall))
            );
            scalar kN = Estar_*LAB;
            //scalar kN =(4.0/3.0)*sqrt(R)*Estar_;

            scalar etaN = alpha_*sqrt(M*kN)*pow025(mag(delta));

            //- Normal velocity
            scalar vnAB(vAB & normal);

            // Normal force
            vector Fn = normal*(kN*pow(delta, b_) - etaN*vnAB);

            // Cohesion force, energy density multiplied by the area of
            // particle-particle overlap
            if (cohesion_)
            {
                Fn -= cohesionEnergyDensity_*0.5*(areaA + areaB)*normal;
            }

            // Apply forcing
            objectA.force()[facei] += Fn;
            objectB.force()[facej] -= Fn;

            vector vtAB = vAB - vnAB*normal;
            vector tangent(vtAB/(mag(vtAB) + small));


            scalar deltaAreaAB = mag(vtAB)*dT;

            areaA += deltaAreaAB;
            areaB -= deltaAreaAB;

            scalar overlapArea = mag(areaA);

            if (overlapArea > vSmall)
            {
                scalar kT = 8.0*sqrt(LAB*delta)*Gstar_;

                scalar etaT = etaN;

                // Tangential force
                vector Ft(Zero);

                if (kT*overlapArea > mu_*mag(Fn))
                {
                    // Tangential force greater than sliding friction,
                    // particle slips
                    Ft = -mu_*mag(Fn)*tangent;
                }
                else
                {
                    Ft = - kT*overlapArea*tangent - etaT*vtAB;
                }

                // Apply forcing
                objectA.force()[facei] += Ft;
                objectB.force()[facej] -= Ft;
            }
        }
    }
}


// ************************************************************************* //
