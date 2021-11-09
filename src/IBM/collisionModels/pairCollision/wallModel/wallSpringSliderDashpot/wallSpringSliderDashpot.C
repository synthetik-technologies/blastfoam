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

#include "wallSpringSliderDashpot.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wallSpringSliderDashpot, 0);
    addToRunTimeSelectionTable
    (
        wallModel,
        wallSpringSliderDashpot,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallSpringSliderDashpot::wallSpringSliderDashpot
(
    const dictionary& dict
)
:
    wallModel(dict, typeName),
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
//     collisionResolutionSteps_
//     (
//         this->coeffDict().template lookup<scalar>("collisionResolutionSteps")
//     ),
{

//     scalar nu = this->coeffDict().template lookup<scalar>("nu");
//
//     scalar E = this->coeffDict().template lookup<scalar>("E");
//
//     scalar pNu = this->owner().constProps().poissonsRatio();
//
//     scalar pE = this->owner().constProps().youngsModulus();
//
//     Estar_ = 1/((1 - sqr(pNu))/pE + (1 - sqr(nu))/E);
//
//     Gstar_ = 1/(2*((2 + pNu - sqr(pNu))/pE + (2 + nu - sqr(nu))/E));
    scalar nu = this->coeffDict().template lookup<scalar>("nu");

    scalar E = this->coeffDict().template lookup<scalar>("E");

    Estar_ = E/(2.0*(1.0 - sqr(nu)));

    scalar G = E/(2.0*(1.0 + nu));

    Gstar_ = G/(2.0*(2.0 - nu));

    cohesion_ = (mag(cohesionEnergyDensity_) > vSmall);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallSpringSliderDashpot::~wallSpringSliderDashpot()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wallSpringSliderDashpot::evaluateWall
(
    const polyPatch& patch,
    immersedBoundaryObject& object,
    const wallCollisionData& map
) const
{
    label facei(map.hitPointIndex());
    if (facei < 0)
    {
        return;
    }

    const scalar dT = object.pMesh().time().deltaTValue();
    const tensor invI = object.invMomentOfInertia();
    const vector& x = object.centre();
    scalar area = 0;
    forAll(map, facej)
    {
        if (map[facej][0] >= 0)
        {
            area += mag(object.Sf()[map[facej][0]]);
        }
    }

    vector hitPoint(map.hitPoint());
    vector normal(x - hitPoint);
    normal = normal/mag(normal);

    const vector& xO = object.faceCentres()[facei];

    scalar delta = (xO - hitPoint) & normal;
    vector F = Zero;

    if (delta > 0)
    {
        // Particles in collision
        vector v(map.v());

        // Effective mass
        vector r = hitPoint - x;
//         scalar M =
//             1.0/(1.0/mass + (normal & (invI & ((r ^ normal) ^ r))));

//         scalar area(object.magSf()[facei]);
        scalar L(area/(mag(r) + vSmall));
        scalar kN = Estar_*L;

        scalar etaN = alpha_;//*sqrt(M*kN)*pow025(delta);

        //- Normal velocity
        scalar vn(v & normal);

        // Normal force
        scalar fn = (kN*pow(delta, b_) - etaN*vn);
        vector Fn = normal*fn;

        // Cohesion force, energy density multiplied by the area of
        // particle-particle overlap
        if (cohesion_)
        {
            Fn -= cohesionEnergyDensity_*area*normal;
        }

        // Apply normal forcing
        F += Fn;

        vector tangent = v - vn*normal;
        tangent = (tangent/(mag(tangent) + small));
        scalar vt = v & tangent;

        scalar overlapArea = mag(area);

//         if (fn < 0)
//         {
//             scalar kT = 8.0*sqrt(L*delta)*Gstar_;
//
//             scalar etaT = etaN;
//
//             // Tangential force
//             vector Ft(Zero);
//
//             if (kT*area > mu_*mag(fn))
//             {
//                 // Tangential force greater than sliding friction,
//                 // particle slips
//                 Ft = -mu_*mag(fn)*tangent;
//             }
//             else
//             {
//                 Ft = - kT*area*tangent - etaT*vt*tangent;
//             }
//
//             // Apply tangential forcing
//             F += Ft;
//         }
    }
//     object.force()[facei] += F;//*map.weights()[facej];
    forAll(map, facej)
    {
        if (map[facej][0] >= 0)
        {
            object.force()[map[facej][0]] += F*map.weights()[facej];
        }
    }
}
//     scalar dT = object.pMesh().time().deltaTValue();
//     scalar mass = object.mass();
//     tensor invI = object.invMomentOfInertia();
//     vector xc = object.centre();
//
//     forAll(map, pairi)
//     {
//         const label facei = map[pairi][0];
//         const label facej = map[pairi][1];
//         if (facei < 0 || facej < 0)
//         {
//             continue;
//         }
//
//         const vector& xO = object.faceCentres()[facei];
//         const vector& xW = map.patch().faceCentres()[facej];
//
//         vector normal(xc - xW);
//         normal /= mag(normal) + vSmall;
//
//         scalar delta = (xO - xW) & normal;
//
//         if (delta > 0)
//         {
//             // Particles in collision
//             vector v(object.velocity(xW));
//
//             // Effective mass
//             vector r(xO - xc);
//             scalar M =
//                 1.0/(1.0/mass + (normal & (invI & ((r ^ normal) ^ r))));
//
//             scalar area(object.magSf()[facei]);
//             scalar L(area/(mag(r) + vSmall));
//             scalar kN = Estar_*L;
//
//             scalar etaN = alpha_;//*sqrt(M*kN)*pow025(delta);
//
//             //- Normal velocity
//             scalar vn(v & normal);
//
//             // Normal force
//             scalar fn = (kN*pow(delta, b_) - etaN*vn);
//             vector Fn = normal*fn;
//
//             // Cohesion force, energy density multiplied by the area of
//             // particle-particle overlap
//             if (cohesion_)
//             {
//                 Fn -= cohesionEnergyDensity_*area*normal;
//             }
//
//             // Apply normal forcing
//             object.force()[facei] += Fn;
//
//             vector tangent = v - vn*normal;
//             tangent = (tangent/(mag(tangent) + small));
//             scalar vt = v & tangent;
//
//             scalar overlapArea = mag(area);
//
//             if (fn < 0)
//             {
//                 scalar kT = 8.0*sqrt(L*delta)*Gstar_;
//
//                 scalar etaT = etaN;
//
//                 // Tangential force
//                 vector Ft(Zero);
//
//                 if (kT*area > mu_*mag(fn))
//                 {
//                     // Tangential force greater than sliding friction,
//                     // particle slips
//                     Ft = -mu_*mag(fn)*tangent;
//                 }
//                 else
//                 {
//                     Ft = - kT*area*tangent - etaT*vt*tangent;
//                 }
//
//                 // Apply tangential forcing
//                 object.force()[facei] += Ft;
//             }
//         }
//     }
// }
// //     scalar kN = (4.0/3.0)*sqrt(pREff)*Estar_;
//     forAll(objectMap, i)
//     {
//         label facei = objectMap[i][0];
//         label patchFacei = objectMap[i][1];
//         const vector& p = object.faceCentres()[facei];
//         const vector& fc = patch.faceCentres()[patchFacei];
//
//         vector r_PW = p - fc;
//
//         vector U_PW = object.velocity(p);// - wallVelocity
//
//         scalar r_PW_mag = mag(r_PW);
//
//         scalar normalOverlapMag = max(-r_PW_mag, 0.0); //(d - r_PW_mag)
//
//         vector rHat_PW = r_PW/(r_PW_mag + vSmall);
//
//         scalar etaN = alpha_*sqrt(object.mass())*pow025(normalOverlapMag);
// //         scalar etaN = alpha_*sqrt(object.mass()*kN)*pow025(normalOverlapMag);
//
//         vector fN_PW =
//             rHat_PW*(- etaN*(U_PW & rHat_PW));
// //         vector fN_PW =
// //             rHat_PW
// //            *(kN*pow(normalOverlapMag, b_) - etaN*(U_PW & rHat_PW));
//
//         // Cohesion force, energy density multiplied by the area of wall/particle
//         // overlap
// //         if (cohesion)
// //         {
// //             fN_PW +=
// //             -cohesionEnergyDensity_
// //             *mathematical::pi*(sqr(pREff) - sqr(r_PW_mag))
// //             *rHat_PW;
// //         }
//
//         object.force()[facei] += fN_PW;
//
//         vector USlip_PW =
//             U_PW - (U_PW & rHat_PW)*rHat_PW
//         + (object.omega() ^ (-rHat_PW)); //(object.omega() ^ (pREff*-rHat_PW))
//
//         scalar deltaT = object.pMesh().time().deltaTValue();
//
// //         vector& tangentialOverlap_PW =
// //             p.collisionRecords().matchWallRecord(-r_PW, pREff).collisionData();
// //
// //         tangentialOverlap_PW += USlip_PW*deltaT;
// //
// //         scalar tangentialOverlapMag = mag(tangentialOverlap_PW);
//
// //         if (tangentialOverlapMag > vSmall)
//         {
// //             scalar kT = 8.0*sqrt(pREff*normalOverlapMag)*Gstar_;
//
//             scalar etaT = etaN;
//
//             // Tangential force
//             vector fT_PW;
//
// //             if (kT*tangentialOverlapMag > mu_*mag(fN_PW))
//             {
//                 // Tangential force greater than sliding friction,
//                 // particle slips
//
//                 fT_PW = -mu_*mag(fN_PW)*USlip_PW/mag(USlip_PW);
//
// //                 tangentialOverlap_PW = Zero;
//             }
// //             else
// //             {
// //                 fT_PW = - kT*tangentialOverlap_PW - etaT*USlip_PW;
// //             }
//
//             object.force()[facei] += fT_PW;
//
// //             p.torque() += (pREff*-rHat_PW) ^ fT_PW;
//         }
//     }
// }


void Foam::wallSpringSliderDashpot::evaluateExternalForce
(
    const polyPatch& patch,
    immersedBoundaryObject& object,
    const wallCollisionData& map
) const
{}

// ************************************************************************* //
