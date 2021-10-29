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

#include "wallLocalSpringSliderDashpot.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wallLocalSpringSliderDashpot, 0);
    addToRunTimeSelectionTable
    (
        wallModel,
        wallLocalSpringSliderDashpot,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallLocalSpringSliderDashpot::wallLocalSpringSliderDashpot
(
    const dictionary& dict
)
:
    wallModel(dict, typeName),
    Estar_(),
    Gstar_(),
    alpha_(),
    b_(),
    mu_(),
    cohesionEnergyDensity_(),
    cohesion_(),
    patchMap_(),
    maxEstarIndex_(-1),
    collisionResolutionSteps_
    (
        this->coeffDict().template lookup<scalar>("collisionResolutionSteps")
    ),
    volumeFactor_(1.0),
    useEquivalentSize_(Switch(this->coeffDict().lookup("useEquivalentSize")))
{
    if (useEquivalentSize_)
    {
        volumeFactor_ =
            this->coeffDict().template lookup<scalar>("volumeFactor");
    }

//     scalar pNu = this->owner().constProps().poissonsRatio();
//
//     scalar pE = this->owner().constProps().youngsModulus();
//
//     const polyMesh& mesh = cloud.mesh();
//
//     const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
//
//     patchMap_.setSize(bMesh.size(), -1);
//
//     DynamicList<label> wallPatchIndices;
//
//     forAll(bMesh, patchi)
//     {
//         if (isA<wallPolyPatch>(bMesh[patchi]))
//         {
//             wallPatchIndices.append(bMesh[patchi].index());
//         }
//     }
//
//     label nWallPatches = wallPatchIndices.size();
//
//     Estar_.setSize(nWallPatches);
//     Gstar_.setSize(nWallPatches);
//     alpha_.setSize(nWallPatches);
//     b_.setSize(nWallPatches);
//     mu_.setSize(nWallPatches);
//     cohesionEnergyDensity_.setSize(nWallPatches);
//     cohesion_.setSize(nWallPatches);
//
//     scalar maxEstar = -great;
//
//     forAll(wallPatchIndices, wPI)
//     {
//         const dictionary& patchCoeffDict
//         (
//             this->coeffDict().subDict(bMesh[wallPatchIndices[wPI]].name())
//         );
//
//         patchMap_[wallPatchIndices[wPI]] = wPI;
//
//         scalar nu = patchCoeffDict.template lookup<scalar>("poissonsRatio");
//
//         scalar E = patchCoeffDict.template lookup<scalar>("youngsModulus");
//
//         Estar_[wPI] = 1/((1 - sqr(pNu))/pE + (1 - sqr(nu))/E);
//
//         Gstar_[wPI] = 1/(2*((2 + pNu - sqr(pNu))/pE + (2 + nu - sqr(nu))/E));
//
//         alpha_[wPI] = patchCoeffDict.template lookup<scalar>("alpha");
//
//         b_[wPI] = patchCoeffDict.template lookup<scalar>("b");
//
//         mu_[wPI] = patchCoeffDict.template lookup<scalar>("mu");
//
//         cohesionEnergyDensity_[wPI] =
//             patchCoeffDict.lookup<scalar>("cohesionEnergyDensity");
//
//         cohesion_[wPI] = (mag(cohesionEnergyDensity_[wPI]) > vSmall);
//
//         if (Estar_[wPI] > maxEstar)
//         {
//             maxEstarIndex_ = wPI;
//
//             maxEstar = Estar_[wPI];
//         }
//     }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallLocalSpringSliderDashpot::~wallLocalSpringSliderDashpot()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wallLocalSpringSliderDashpot::evaluateWall
(
    const polyPatch& patch,
    immersedBoundaryObject& object,
    const wallCollisionData& map
) const
{
//     // wall patch index
//     label wPI = patchMap_[data.patchIndex()];
//
//     // data for this patch
//     scalar Estar = Estar_[wPI];
//     scalar Gstar = Gstar_[wPI];
//     scalar alpha = alpha_[wPI];
//     scalar b = b_[wPI];
//     scalar mu = mu_[wPI];
//     scalar cohesionEnergyDensity = cohesionEnergyDensity_[wPI];
//     cohesion = cohesion && cohesion_[wPI];
//
//     vector r_PW = p.position() - site;
//
//     vector U_PW = object.U() - data.wallData();
//
//     scalar r_PW_mag = mag(r_PW);
//
//     scalar normalOverlapMag = max(pREff - r_PW_mag, 0.0);
//
//     vector rHat_PW = r_PW/(r_PW_mag + vSmall);
//
//     scalar kN = (4.0/3.0)*sqrt(pREff)*Estar;
//
//     scalar etaN = alpha*sqrt(object.mass()*kN)*pow025(normalOverlapMag);
//
//     vector fN_PW =
//         rHat_PW
//        *(kN*pow(normalOverlapMag, b) - etaN*(U_PW & rHat_PW));
//
//     // Cohesion force, energy density multiplied by the area of wall/particle
//     // overlap
//     if (cohesion)
//     {
//         fN_PW +=
//            -cohesionEnergyDensity
//            *mathematical::pi*(sqr(pREff) - sqr(r_PW_mag))
//            *rHat_PW;
//     }
//
//     p.f() += fN_PW;
//
//     vector USlip_PW =
//         U_PW - (U_PW & rHat_PW)*rHat_PW
//       + (p.omega() ^ (pREff*-rHat_PW));
//
//     scalar deltaT = this->owner().mesh().time().deltaTValue();
//
//     vector& tangentialOverlap_PW =
//         p.collisionRecords().matchWallRecord(-r_PW, pREff).collisionData();
//
//     tangentialOverlap_PW += USlip_PW*deltaT;
//
//     scalar tangentialOverlapMag = mag(tangentialOverlap_PW);
//
//     if (tangentialOverlapMag > vSmall)
//     {
//         scalar kT = 8.0*sqrt(pREff*normalOverlapMag)*Gstar;
//
//         scalar etaT = etaN;
//
//         // Tangential force
//         vector fT_PW;
//
//         if (kT*tangentialOverlapMag > mu*mag(fN_PW))
//         {
//             // Tangential force greater than sliding friction,
//             // particle slips
//
//             fT_PW = -mu*mag(fN_PW)*USlip_PW/mag(USlip_PW);
//
//             tangentialOverlap_PW = Zero;
//         }
//         else
//         {
//             fT_PW =
//                 -kT*tangentialOverlapMag
//                *tangentialOverlap_PW/tangentialOverlapMag
//               - etaT*USlip_PW;
//         }
//
//         p.f() += fT_PW;
//
//         p.torque() += (pREff*-rHat_PW) ^ fT_PW;
//     }
}


void Foam::wallLocalSpringSliderDashpot::evaluateExternalForce
(
    const polyPatch& patch,
    immersedBoundaryObject& object,
    const wallCollisionData& map
) const
{}


// ************************************************************************* //
