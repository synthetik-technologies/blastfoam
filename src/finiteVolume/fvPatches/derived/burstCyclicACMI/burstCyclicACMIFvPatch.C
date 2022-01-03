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

#include "burstCyclicACMIFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(burstCyclicACMIFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, burstCyclicACMIFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::burstCyclicACMIFvPatch::delta() const
{
    return
        intact_*fvPatch::delta()
      + (1.0 - intact_)*cyclicACMIFvPatch::delta();
}


void Foam::burstCyclicACMIFvPatch::autoMap
(
    const fvPatchFieldMapper& m
)
{
    if (!burstCyclicACMIPolyPatch_.needMap())
    {
        return;
    }
    m(intact_, intact_);
    intact_ = pos(intact_ - 0.5);
}


void Foam::burstCyclicACMIFvPatch::rmap
(
    const fvPatch& ptf,
    const labelList& addr
)
{
    if (!burstCyclicACMIPolyPatch_.needMap())
    {
        return;
    }
    const burstCyclicACMIFvPatch& bp = refCast<const burstCyclicACMIFvPatch>(ptf);
    intact_.rmap(bp.intact_, addr);
    intact_ = pos(intact_ - 0.5);
}


void Foam::burstCyclicACMIFvPatch::update
(
    const scalarField& p,
    const scalarField& impulse
)
{
    if (burstCyclicACMIPolyPatch_.partialBurst())
    {
        scalarField op(intact_.size(), -great);
        scalarField imp(intact_.size(), -great);
        if (burstCyclicACMIPolyPatch_.usePressure())
        {
            op = p - burstCyclicACMIPolyPatch_.pBurst();
        }
        if (burstCyclicACMIPolyPatch_.useImpulse())
        {
            imp = impulse - burstCyclicACMIPolyPatch_.impulseBurst();
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
        // Patch has already burstCyclicACMI
        if (intact_[0] == 0)
        {
            return;
        }

        if (burstCyclicACMIPolyPatch_.usePressure())
        {
            scalar maxP(max(p));
            intact_ = maxP > burstCyclicACMIPolyPatch_.pBurst() ? 0.0 : 1.0;
        }
        if (burstCyclicACMIPolyPatch_.useImpulse())
        {
            scalar maxImpulse(max(impulse));
            intact_ = maxImpulse > burstCyclicACMIPolyPatch_.impulseBurst() ? 0.0 : 1.0;
        }
    }
}

// ************************************************************************* //
