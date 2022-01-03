/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "burstFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(burstFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, burstFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::burstFvPatch::autoMap
(
    const fvPatchFieldMapper& m
)
{
    if (!burstPolyPatch_.needMap())
    {
        return;
    }
    m(intact_, intact_);
    intact_ = pos(intact_ - 0.5);

    if (burstPolyPatch_.usePressure())
    {
        m(pRef_, pRef_);
    }
}


void Foam::burstFvPatch::rmap
(
    const fvPatch& ptf,
    const labelList& addr
)
{
    if (!burstPolyPatch_.needMap())
    {
        return;
    }
    const burstFvPatch& bp = refCast<const burstFvPatch>(ptf);
    intact_.rmap(bp.intact_, addr);
    intact_ = pos(intact_ - 0.5);

    if (burstPolyPatch_.usePressure())
    {
        pRef_.rmap(bp.pRef_, addr);
    }
}


void Foam::burstFvPatch::update
(
    const scalarField& p,
    const scalarField& impulse
)
{
    if (burstPolyPatch_.partialBurst())
    {
        scalarField op(intact_.size(), -great);
        scalarField imp(intact_.size(), -great);
        if (burstPolyPatch_.usePressure())
        {
            op = p - burstPolyPatch_.pBurst() - pRef_;
        }
        if (burstPolyPatch_.useImpulse())
        {
            imp = impulse - burstPolyPatch_.impulseBurst();
        }
        forAll(intact_, facei)
        {
            if (intact_[facei])
            {
                intact_[facei] = op[facei] < 0 && imp[facei] < 0;
            }
        }
    }
    else if (!p.size())
    {
        return;
    }
    else
    {
        // Patch has already burst
        if (intact_[0] == 0)
        {
            return;
        }

        if (burstPolyPatch_.usePressure())
        {
            scalar maxP(max(p));
            intact_ = maxP > burstPolyPatch_.pBurst() ? 0.0 : 1.0;
        }
        if (burstPolyPatch_.useImpulse())
        {
            scalar maxImpulse(max(impulse));
            intact_ = maxImpulse > burstPolyPatch_.impulseBurst() ? 0.0 : 1.0;
        }
    }
}

// ************************************************************************* //
