/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024
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

#include "LimitedScheme.H"
#include "Limited01.H"
#include "logarithmic.H"

#include "ReconstructionScheme.H"
#include "MUSCLReconstructionScheme.H"
#include "LinearMUSCLReconstructionScheme.H"
#include "QuadraticMUSCLReconstructionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeLimitedSurfaceInterpolationScheme(logarithmic, logarithmicLimiter)
    makeLimitedVSurfaceInterpolationScheme(logarithmicV,logarithmicLimiter)

    makeLLimitedSurfaceInterpolationTypeScheme
    (
        limitedLogarithmic,
        LimitedLimiter,
        logarithmicLimiter,
        NVDTVD,
        magSqr,
        scalar
    )

    makeLLimitedSurfaceInterpolationTypeScheme
    (
        logarithmic01,
        Limited01Limiter,
        logarithmicLimiter,
        NVDTVD,
        magSqr,
        scalar
    )

    declareReconstructionSchemes(LinearMUSCL)
    makeMUSCLReconstruction(Linear, "linearMUSCL", logarithmic, logarithmicLimiter)
    makeLMUSCLReconstruction
    (
        Linear,
        "linearMUSCL",
        logarithmic01,
        Limited01Limiter,
        logarithmicLimiter,
        NVDTVD,
        magSqr,
        scalar
    )

    declareReconstructionSchemes(QuadraticMUSCL)
    makeMUSCLReconstruction(Quadratic, "quadraticMUSCL", logarithmic, logarithmicLimiter)
    makeLMUSCLReconstruction
    (
        Quadratic,
        "quadraticMUSCL",
        logarithmic01,
        Limited01Limiter,
        logarithmicLimiter,
        NVDTVD,
        magSqr,
        scalar
    )
}

// ************************************************************************* //
