/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "adaptiveTypesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeAdaptiveType
(
    1,
    mag(f),
    f,
    mag(f1 - f2)/stabilise(mag(f1), small),
    scalar
);
makeAdaptiveType
(
    pTraits<vector>::nComponents,
    mag(f),
    f[i],
    mag(f1 - f2)/stabilise(mag(f1), small),
    vector
);
makeAdaptiveType
(
    pTraits<symmTensor>::nComponents,
    mag(f),
    f[i],
    mag(f1 - f2)/stabilise(mag(f1), small),
    symmTensor
);
makeAdaptiveType
(
    pTraits<sphericalTensor>::nComponents,
    mag(f),
    f[i],
    mag(f1 - f2)/stabilise(mag(f1), small),
    sphericalTensor
);
makeAdaptiveType
(
    pTraits<tensor>::nComponents,
    mag(f),
    f[i],
    mag(f1 - f2)/stabilise(mag(f1), small),
    tensor
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}

// ************************************************************************* //
