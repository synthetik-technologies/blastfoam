/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "LimitedScheme.H"
#include "Limited01.H"
#include "overbee.H"

#include "ReconstructionScheme.H"
#include "MUSCLReconstructionScheme.H"
#include "LinearMUSCLReconstructionScheme.H"
#include "QuadraticMUSCLReconstructionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeLimitedSurfaceInterpolationScheme(overbee, overbeeLimiter)
    makeLimitedVSurfaceInterpolationScheme(overbeeV, overbeeLimiter)

    makeLLimitedSurfaceInterpolationTypeScheme
    (
        limitedoverbee,
        LimitedLimiter,
        overbeeLimiter,
        NVDTVD,
        magSqr,
        scalar
    )

    makeLLimitedSurfaceInterpolationTypeScheme
    (
        overbee01,
        Limited01Limiter,
        overbeeLimiter,
        NVDTVD,
        magSqr,
        scalar
    );

    declareReconstructionSchemes(LinearMUSCL)
    makeMUSCLReconstruction(Linear, "linearMUSCL", overbee, overbeeLimiter)
    makeLMUSCLReconstruction
    (
        Linear,
        "linearMUSCL",
        overbee01,
        Limited01Limiter,
        overbeeLimiter,
        NVDTVD,
        magSqr,
        scalar
    )

    declareReconstructionSchemes(QuadraticMUSCL)
    makeMUSCLReconstruction(Quadratic, "quadraticMUSCL", overbee, overbeeLimiter)
    makeLMUSCLReconstruction
    (
        Quadratic,
        "quadraticMUSCL",
        overbee01,
        Limited01Limiter,
        overbeeLimiter,
        NVDTVD,
        magSqr,
        scalar
    )
}

// ************************************************************************* //
