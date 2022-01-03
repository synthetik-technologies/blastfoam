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

#include "burstCyclicAMIFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(burstCyclicAMIFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, burstCyclicAMIFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::burstCyclicAMIFvPatch::delta() const
{
   return
        intact_*fvPatch::delta()
      + (1.0 - intact_)*cyclicAMIFvPatch::delta();
}


void Foam::burstCyclicAMIFvPatch::autoMap
(
    const fvPatchFieldMapper& m
)
{
    if (!burstCyclicAMIPolyPatch_.needMap())
    {
        return;
    }
    m(intact_, intact_);
    intact_ = pos(intact_ - 0.5);
}


void Foam::burstCyclicAMIFvPatch::rmap
(
    const fvPatch& ptf,
    const labelList& addr
)
{
    if (!burstCyclicAMIPolyPatch_.needMap())
    {
        return;
    }
    const burstCyclicAMIFvPatch& bp = refCast<const burstCyclicAMIFvPatch>(ptf);
    intact_.rmap(bp.intact_, addr);
    intact_ = pos(intact_ - 0.5);
}


void Foam::burstCyclicAMIFvPatch::update
(
    const scalarField& p,
    const scalarField& impulse
)
{
    if (burstCyclicAMIPolyPatch_.partialBurst())
    {
        scalarField op(intact_.size(), -great);
        scalarField imp(intact_.size(), -great);
        if (burstCyclicAMIPolyPatch_.usePressure())
        {
            op = p - burstCyclicAMIPolyPatch_.pBurst();
        }
        if (burstCyclicAMIPolyPatch_.useImpulse())
        {
            imp = impulse - burstCyclicAMIPolyPatch_.impulseBurst();
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
        // Patch has already burstCyclicAMI
        if (intact_[0] == 0)
        {
            return;
        }

        if (burstCyclicAMIPolyPatch_.usePressure())
        {
            scalar maxP(max(p));
            intact_ = maxP > burstCyclicAMIPolyPatch_.pBurst() ? 0.0 : 1.0;
        }
        if (burstCyclicAMIPolyPatch_.useImpulse())
        {
            scalar maxImpulse(max(impulse));
            intact_ = maxImpulse > burstCyclicAMIPolyPatch_.impulseBurst() ? 0.0 : 1.0;
        }
    }
}

// ************************************************************************* //
