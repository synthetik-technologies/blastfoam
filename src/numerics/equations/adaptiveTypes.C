/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
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
    Foam::mag(f),
    f,
    Foam::mag(f1 - f2),
    scalar
);
makeAdaptiveType
(
    pTraits<vector>::nComponents,
    Foam::mag(f),
    f[i],
    Foam::mag(f1 - f2),
    vector
);
makeAdaptiveType
(
    pTraits<symmTensor>::nComponents,
    Foam::mag(f),
    f[i],
    Foam::mag(f1 - f2),
    symmTensor
);
makeAdaptiveType
(
    pTraits<sphericalTensor>::nComponents,
    Foam::mag(f),
    f[i],
    Foam::mag(f1 - f2),
    sphericalTensor
);
makeAdaptiveType
(
    pTraits<tensor>::nComponents,
    Foam::mag(f),
    f[i],
    Foam::mag(f1 - f2),
    tensor
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}

// ************************************************************************* //
