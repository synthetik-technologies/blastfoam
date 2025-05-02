/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dynamicRelaxation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicRelaxation, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicRelaxation::dynamicRelaxation
(
    const dictionary& dict
)
:
    relax_(false),
    start_(-great),
    nSteps_(1),
    curIndex_(-1),
    Ks_(0),
    scale_(0.0),
    useProbe_(false),
    probe_(Zero),
    curScale_(1.0)
{
    read(dict);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dynamicRelaxation::read(const dictionary& dict)
{
    dict.readIfPresent("dynamicRelaxation", relax_);
    if (relax_)
    {
        dict.readIfPresent("startRelaxation", start_);
        dict.readIfPresent("nRelaxationSteps", nSteps_);
        Ks_.setCapacity(nSteps_);
        dict.readIfPresent("relaxationScale", scale_);

        dict.readIfPresent("useProbe", useProbe_);
        if (useProbe_)
        {
            probe_ = dict.lookup<vector>("probeLocation");
        }
    }
}


bool Foam::dynamicRelaxation::relax
(
    volVectorField& U,
    const volScalarField& rho
)
{
    if (!relax_ || U.time().value() < start_)
    {
        curScale_ = 1.0;
        return false;
    }
    if (curIndex_ != U.time().timeIndex())
    {
        if (Ks_.size() < nSteps_)
        {
            Ks_.append(0.0);
        }
        for (label i = Ks_.size() - 1; i > 0; i--)
        {
            Ks_[i] = Ks_[i-1];
        }
        curIndex_ = U.time().timeIndex();
    }

    scalar K = 0.0;
    if (useProbe_)
    {
        label celli = U.mesh().findNearestCell(probe_);
        scalar dist = magSqr(U.mesh().C()[celli] - probe_);
        scalar minDist = returnReduce(dist, minOp<scalar>());
        if (dist == minDist)
        {
            K = 0.5*magSqr(U[celli])*rho[celli];
        }
        reduce(K, maxOp<scalar>());

        if (returnReduce(celli, maxOp<label>()) == -1)
        {
            WarningInFunction
                << "Could not find cell with location " << probe_ << nl
                << "Using total kinetic energy" << endl;
        }
    }
    else
    {
        K =
            gSum
            (
                0.5*magSqr(U.internalField())
               *rho.primitiveField()
               *U.mesh().V().primitiveField()
            );
    }
    Ks_[0] = K;

    if (Ks_.size() < nSteps_)
    {
        curScale_ = 1.0;
        return false;
    }

    bool limit = true;
    for (label i = 1; i < Ks_.size(); i++)
    {
        if (Ks_[i] < K)
        {
            limit = false;
            break;
        }
    }
    if (limit)
    {
        curScale_ = scale_;
        U *= curScale_;
        Ks_.clear();
        return true;
    }

    curScale_ = 1.0;
    return false;
}

// ************************************************************************* //
