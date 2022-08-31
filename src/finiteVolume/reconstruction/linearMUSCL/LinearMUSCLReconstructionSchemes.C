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
#include "LinearMUSCLReconstructionScheme.H"


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
defineReconstructionTable(Linear, LinearMUSCL, scalar);
defineReconstructionTable(Linear, LinearMUSCL, vector);
defineReconstructionTable(Linear, LinearMUSCL, symmTensor);
defineReconstructionTable(Linear, LinearMUSCL, sphericalTensor);
defineReconstructionTable(Linear, LinearMUSCL, tensor);

// Define limiters
makeMUSCLReconstruction(Linear, "linearMUSCL", Gamma, GammaLimiter);
makeLMUSCLReconstruction
(
    Linear,
    "linearMUSCL",
    limitedGamma,
    LimitedLimiter,
    GammaLimiter,
    NVDTVD,
    magSqr,
    scalar
);
makeLMUSCLReconstruction
(
    Linear,
    "linearMUSCL",
    Gamma01,
    Limited01Limiter,
    GammaLimiter,
    NVDTVD,
    magSqr,
    scalar
);

makeMUSCLReconstruction(Linear, "linearMUSCL", MUSCL, MUSCLLimiter);
makeLMUSCLReconstruction
(
    Linear,
    "linearMUSCL",
    limitedMUSCL,
    LimitedLimiter,
    MUSCLLimiter,
    NVDTVD,
    magSqr,
    scalar
);
makeLMUSCLReconstruction
(
    Linear,
    "linearMUSCL",
    MUSCL01,
    Limited01Limiter,
    MUSCLLimiter,
    NVDTVD,
    magSqr,
    scalar
);

makeMUSCLReconstruction(Linear, "linearMUSCL", Minmod, MinmodLimiter);

makeMUSCLReconstruction(Linear, "linearMUSCL", OSPRE, OSPRELimiter);

makeMUSCLReconstruction(Linear, "linearMUSCL", QUICK, QUICKLimiter);

makeMUSCLReconstruction(Linear, "linearMUSCL", SFCD, SFCDLimiter);

makeMUSCLReconstruction(Linear, "linearMUSCL", SuperBee, SuperBeeLimiter);

makeMUSCLReconstruction(Linear, "linearMUSCL", UMIST, UMISTLimiter);

makeMUSCLReconstruction(Linear, "linearMUSCL", limitedCubic, limitedCubicLimiter);
makeLMUSCLReconstruction
(
    Linear,
    "linearMUSCL",
    limitedLimitedCubic,
    LimitedLimiter,
    limitedCubicLimiter,
    NVDTVD,
    magSqr,
    scalar
);
makeLMUSCLReconstruction
(
    Linear,
    "linearMUSCL",
    limitedCubic01,
    Limited01Limiter,
    limitedCubicLimiter,
    NVDTVD,
    magSqr,
    scalar
);

makeMUSCLReconstruction(Linear, "linearMUSCL", limitedLinear, limitedLinearLimiter);
makeLMUSCLReconstruction
(
    Linear,
    "linearMUSCL",
    limitedLimitedLinear,
    LimitedLimiter,
    limitedLinearLimiter,
    NVDTVD,
    magSqr,
    scalar
);
makeLMUSCLReconstruction
(
    Linear,
    "linearMUSCL",
    limitedLinear01,
    Limited01Limiter,
    limitedLinearLimiter,
    NVDTVD,
    magSqr,
    scalar
);

makeMUSCLReconstruction(Linear, "linearMUSCL", vanAlbada, vanAlbadaLimiter);

makeMUSCLReconstruction(Linear, "linearMUSCL", vanLeer, vanLeerLimiter);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
