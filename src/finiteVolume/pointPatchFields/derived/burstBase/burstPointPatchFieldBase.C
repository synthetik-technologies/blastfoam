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

#include "burstPointPatchFieldBase.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(burstPointPatchFieldBase, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::burstPointPatchFieldBase::burstPointPatchFieldBase
(
    const polyPatch& pp
)
:
    burstBase_
    (
        const_cast<burstPolyPatchBase&>
        (
            dynamicCast<const burstPolyPatchBase>(pp)
        )
    )
{}


Foam::burstPointPatchFieldBase::burstPointPatchFieldBase
(
    const polyPatch& pp,
    const burstPointPatchFieldBase& bpf
)
:
    burstBase_
    (
        const_cast<burstPolyPatchBase&>
        (
            dynamicCast<const burstPolyPatchBase>(pp)
        )
    )
{}


Foam::burstPointPatchFieldBase::burstPointPatchFieldBase
(
    const polyPatch& pp,
    const burstPointPatchFieldBase& bpf,
    const pointPatchFieldMapper& mapper
)
:
    burstBase_
    (
        const_cast<burstPolyPatchBase&>
        (
            dynamicCast<const burstPolyPatchBase>(pp)
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burstPointPatchFieldBase::~burstPointPatchFieldBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
