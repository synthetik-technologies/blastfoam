/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020
     \\/     M anipulation  | Synthetik Applied Technology
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

\*---------------------------------------------------------------------------*/

#include "ReconstructionScheme.H"
#include "MUSCLReconstructionScheme.H"
#include "QuadraticMUSCLReconstructionScheme.H"


#include "Gamma.H"
#include "Limited01.H"
#include "MUSCL.H"
#include "Minmod.H"
#include "OSPRE.H"
#include "QUICK.H"
#include "SFCD.H"
#include "SuperBee.H"
#include "UMIST.H"
#include "limitedCubic.H"
#include "limitedLinear.H"
#include "vanLeer.H"
#include "vanAlbada.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Define runtime selection tables
defineReconstructionTable(Quadratic, QuadraticMUSCL, scalar)
defineReconstructionTable(Quadratic, QuadraticMUSCL, vector)
defineReconstructionTable(Quadratic, QuadraticMUSCL, symmTensor)
defineReconstructionTable(Quadratic, QuadraticMUSCL, sphericalTensor)
defineReconstructionTable(Quadratic, QuadraticMUSCL, tensor)

// Define limiters
makeMUSCLReconstruction(Quadratic, "quadraticMUSCL", Gamma, GammaLimiter)
makeLMUSCLReconstruction
(
    Quadratic,
    "quadraticMUSCL",
    limitedGamma,
    LimitedLimiter,
    GammaLimiter,
    NVDTVD,
    magSqr,
    scalar
)
makeLMUSCLReconstruction
(
    Quadratic,
    "quadraticMUSCL",
    Gamma01,
    Limited01Limiter,
    GammaLimiter,
    NVDTVD,
    magSqr,
    scalar
)

makeMUSCLReconstruction(Quadratic, "quadraticMUSCL", MUSCL, MUSCLLimiter)
makeLMUSCLReconstruction
(
    Quadratic, "quadraticMUSCL",
    limitedMUSCL,
    LimitedLimiter,
    MUSCLLimiter,
    NVDTVD,
    magSqr,
    scalar
)
makeLMUSCLReconstruction
(
    Quadratic, "quadraticMUSCL",
    MUSCL01,
    Limited01Limiter,
    MUSCLLimiter,
    NVDTVD,
    magSqr,
    scalar
)

makeMUSCLReconstruction(Quadratic, "quadraticMUSCL", Minmod, MinmodLimiter)

makeMUSCLReconstruction(Quadratic, "quadraticMUSCL", OSPRE, OSPRELimiter)

makeMUSCLReconstruction(Quadratic, "quadraticMUSCL", QUICK, QUICKLimiter)

makeMUSCLReconstruction(Quadratic, "quadraticMUSCL", SFCD, SFCDLimiter)

makeMUSCLReconstruction(Quadratic, "quadraticMUSCL", SuperBee, SuperBeeLimiter)

makeMUSCLReconstruction(Quadratic, "quadraticMUSCL", UMIST, UMISTLimiter)

makeMUSCLReconstruction(Quadratic, "quadraticMUSCL", limitedCubic, limitedCubicLimiter)
makeLMUSCLReconstruction
(
    Quadratic, "quadraticMUSCL",
    limitedLimitedCubic,
    LimitedLimiter,
    limitedCubicLimiter,
    NVDTVD,
    magSqr,
    scalar
)
makeLMUSCLReconstruction
(
    Quadratic, "quadraticMUSCL",
    limitedCubic01,
    Limited01Limiter,
    limitedCubicLimiter,
    NVDTVD,
    magSqr,
    scalar
)

makeMUSCLReconstruction(Quadratic, "quadraticMUSCL", limitedLinear, limitedLinearLimiter)
makeLMUSCLReconstruction
(
    Quadratic, "quadraticMUSCL",
    limitedLimitedLinear,
    LimitedLimiter,
    limitedLinearLimiter,
    NVDTVD,
    magSqr,
    scalar
)
makeLMUSCLReconstruction
(
    Quadratic, "quadraticMUSCL",
    limitedLinear01,
    Limited01Limiter,
    limitedLinearLimiter,
    NVDTVD,
    magSqr,
    scalar
)

makeMUSCLReconstruction(Quadratic, "quadraticMUSCL", vanAlbada, vanAlbadaLimiter)

makeMUSCLReconstruction(Quadratic, "quadraticMUSCL", vanLeer, vanLeerLimiter)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
