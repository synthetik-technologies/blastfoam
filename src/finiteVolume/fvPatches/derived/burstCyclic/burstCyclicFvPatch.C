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

#include "burstCyclicFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(burstCyclicFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, burstCyclicFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::burstCyclicFvPatch::makeWeights(scalarField& w) const
{
    const cyclicFvPatch& nbrPatch = neighbFvPatch();

    const vectorField delta(coupledFvPatch::delta());
    const vectorField nbrDelta(nbrPatch.coupledFvPatch::delta());

    const scalarField nfDelta(nf() & delta);
    const scalarField nbrNfDelta(nbrPatch.nf() & nbrDelta);

    forAll(delta, facei)
    {
        if (intact_[facei])
        {
            w[facei] = 1.0;
        }
        else
        {
            const scalar ndoi = nfDelta[facei];
            const scalar ndni = nbrNfDelta[facei];
            const scalar ndi = ndoi + ndni;

            if (ndni/vGreat < ndi)
            {
                w[facei] = ndni/ndi;
            }
            else
            {
                const scalar doi = mag(delta[facei]);
                const scalar dni = mag(nbrDelta[facei]);
                const scalar di = doi + dni;

                w[facei] = dni/di;
            }
        }
    }
}


Foam::tmp<Foam::vectorField> Foam::burstCyclicFvPatch::delta() const
{
    const vectorField D(fvPatch::delta());
    const vectorField patchD(coupledFvPatch::delta());
    const vectorField nbrPatchD(neighbFvPatch().coupledFvPatch::delta());

    tmp<vectorField> tpdv(new vectorField(patchD.size()));
    vectorField& pdv = tpdv.ref();

    // To the transformation if necessary
    if (transform().transforms())
    {
        forAll(patchD, facei)
        {
            if (intact_[facei] > 0.5)
            {
                pdv[facei] = D[facei];
            }
            else
            {
                vector ddi = patchD[facei];
                vector dni = nbrPatchD[facei];

                pdv[facei] = ddi - transform().transform(dni);
            }
        }
    }
    else
    {
        forAll(patchD, facei)
        {
            if (intact_[facei] > 0.5)
            {
                pdv[facei] = D[facei];
            }
            else
            {
                vector ddi = patchD[facei];
                vector dni = nbrPatchD[facei];

                pdv[facei] = ddi - dni;
            }
        }
    }
    return tpdv;
}


void Foam::burstCyclicFvPatch::autoMap
(
    const fvPatchFieldMapper& m
)
{
    if (!burstCyclicPolyPatch_.needMap())
    {
        return;
    }
    m(intact_, intact_);
    intact_ = pos(intact_ - 0.5);
}


void Foam::burstCyclicFvPatch::rmap
(
    const fvPatch& ptf,
    const labelList& addr
)
{
    if (!burstCyclicPolyPatch_.needMap())
    {
        return;
    }
    const burstCyclicFvPatch& bp = refCast<const burstCyclicFvPatch>(ptf);
    intact_.rmap(bp.intact_, addr);
    intact_ = pos(intact_ - 0.5);
}


void Foam::burstCyclicFvPatch::update
(
    const scalarField& p,
    const scalarField& impulse
)
{
    if (burstCyclicPolyPatch_.partialBurst())
    {
        scalarField op(intact_.size(), -great);
        scalarField imp(intact_.size(), -great);
        if (burstCyclicPolyPatch_.usePressure())
        {
            op = p - burstCyclicPolyPatch_.pBurst();
        }
        if (burstCyclicPolyPatch_.useImpulse())
        {
            imp = impulse - burstCyclicPolyPatch_.impulseBurst();
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
        // Patch has already burstCyclic
        if (intact_[0] == 0)
        {
            return;
        }

        if (burstCyclicPolyPatch_.usePressure())
        {
            scalar maxP(max(p));
            intact_ = maxP > burstCyclicPolyPatch_.pBurst() ? 0.0 : 1.0;
        }
        if (burstCyclicPolyPatch_.useImpulse())
        {
            scalar maxImpulse(max(impulse));
            intact_ = maxImpulse > burstCyclicPolyPatch_.impulseBurst() ? 0.0 : 1.0;
        }
    }
}

// ************************************************************************* //
