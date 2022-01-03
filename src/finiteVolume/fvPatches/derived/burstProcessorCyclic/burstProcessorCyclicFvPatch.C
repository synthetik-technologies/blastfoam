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

#include "burstProcessorCyclicFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(burstProcessorCyclicFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, burstProcessorCyclicFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::burstProcessorCyclicFvPatch::makeWeights(scalarField& w) const
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


Foam::tmp<Foam::vectorField> Foam::burstProcessorCyclicFvPatch::delta() const
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


void Foam::burstProcessorCyclicFvPatch::autoMap
(
    const fvPatchFieldMapper& m
)
{
    if (!burstProcessorCyclicPolyPatch_.needMap())
    {
        return;
    }
    m(intact_, intact_);
    intact_ = pos(intact_ - 0.5);
}


void Foam::burstProcessorCyclicFvPatch::rmap
(
    const fvPatch& ptf,
    const labelList& addr
)
{
    if (!burstProcessorCyclicPolyPatch_.needMap())
    {
        return;
    }
    const burstProcessorCyclicFvPatch& bp = refCast<const burstProcessorCyclicFvPatch>(ptf);
    intact_.rmap(bp.intact_, addr);
    intact_ = pos(intact_ - 0.5);
}


void Foam::burstProcessorCyclicFvPatch::update
(
    const scalarField& p,
    const scalarField& impulse
)
{
    if (burstProcessorCyclicPolyPatch_.partialBurst())
    {
        scalarField op(intact_.size(), -great);
        scalarField imp(intact_.size(), -great);
        if (burstProcessorCyclicPolyPatch_.usePressure())
        {
            op = p - burstProcessorCyclicPolyPatch_.pBurst();
        }
        if (burstProcessorCyclicPolyPatch_.useImpulse())
        {
            imp = impulse - burstProcessorCyclicPolyPatch_.impulseBurst();
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
        // Patch has already burstProcessorCyclic
        if (intact_[0] == 0)
        {
            return;
        }

        if (burstProcessorCyclicPolyPatch_.usePressure())
        {
            scalar maxP(max(p));
            intact_ = maxP > burstProcessorCyclicPolyPatch_.pBurst() ? 0.0 : 1.0;
        }
        if (burstProcessorCyclicPolyPatch_.useImpulse())
        {
            scalar maxImpulse(max(impulse));
            intact_ = maxImpulse > burstProcessorCyclicPolyPatch_.impulseBurst() ? 0.0 : 1.0;
        }
    }
}

// ************************************************************************* //
