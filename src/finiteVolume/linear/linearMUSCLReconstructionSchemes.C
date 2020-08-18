/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
18-08-2020 Jeff Heylmun:    | MUSCL reconstruction
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

#include "MUSCLReconstructionScheme.H"
#include "MUSCLReconstruction.H"
#include "linearMUSCLReconstructionScheme.H"


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


// Define types and runtime selection tables
defineMUSCLReconstructionTypes(linear, scalar)
defineMUSCLReconstructionTypes(linear, vector)
defineMUSCLReconstructionTypes(linear, symmTensor)
defineMUSCLReconstructionTypes(linear, sphericalTensor)
defineMUSCLReconstructionTypes(linear, tensor)

// Define limiters
makeMUSCLReconstruction(linear, Gamma, GammaLimiter)
makeLMUSCLReconstruction
(
    linear,
    limitedGamma,
    LimitedLimiter,
    GammaLimiter,
    NVDTVD,
    magSqr,
    scalar
)
makeLMUSCLReconstruction
(
    linear,
    Gamma01,
    Limited01Limiter,
    GammaLimiter,
    NVDTVD,
    magSqr,
    scalar
)

makeMUSCLReconstruction(linear, MUSCL, MUSCLLimiter)
makeLMUSCLReconstruction
(
    linear,
    limitedMUSCL,
    LimitedLimiter,
    MUSCLLimiter,
    NVDTVD,
    magSqr,
    scalar
)
makeLMUSCLReconstruction
(
    linear,
    MUSCL01,
    Limited01Limiter,
    MUSCLLimiter,
    NVDTVD,
    magSqr,
    scalar
)

makeMUSCLReconstruction(linear, Minmod, MinmodLimiter)

makeMUSCLReconstruction(linear, OSPRE, OSPRELimiter)

makeMUSCLReconstruction(linear, QUICK, QUICKLimiter)

makeMUSCLReconstruction(linear, SFCD, SFCDLimiter)

makeMUSCLReconstruction(linear, SuperBee, SuperBeeLimiter)

makeMUSCLReconstruction(linear, UMIST, UMISTLimiter)

makeMUSCLReconstruction(linear, limitedCubic, limitedCubicLimiter)
makeLMUSCLReconstruction
(
    linear,
    limitedLimitedCubic,
    LimitedLimiter,
    limitedCubicLimiter,
    NVDTVD,
    magSqr,
    scalar
)
makeLMUSCLReconstruction
(
    linear,
    limitedCubic01,
    Limited01Limiter,
    limitedCubicLimiter,
    NVDTVD,
    magSqr,
    scalar
)

makeMUSCLReconstruction(linear, limitedLinear, limitedLinearLimiter)
makeLMUSCLReconstruction
(
    linear,
    limitedLimitedLinear,
    LimitedLimiter,
    limitedLinearLimiter,
    NVDTVD,
    magSqr,
    scalar
)
makeLMUSCLReconstruction
(
    linear,
    limitedLinear01,
    Limited01Limiter,
    limitedLinearLimiter,
    NVDTVD,
    magSqr,
    scalar
)

makeMUSCLReconstruction(linear, vanAlbada, vanAlbadaLimiter)

makeMUSCLReconstruction(linear, vanLeer, vanLeerLimiter)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
