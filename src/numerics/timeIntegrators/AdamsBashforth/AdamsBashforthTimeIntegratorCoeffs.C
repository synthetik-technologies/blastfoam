/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020-2023
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "AdamsBashforthTimeIntegratorCoeffs.H"
#include "Field.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace timeIntegrators
{
    defineTypeNameAndDebug(AdamsBashforth, 0);
    addToRunTimeSelectionTable(timeIntegratorCoeffs, AdamsBashforth, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeIntegrators::AdamsBashforth::AdamsBashforth(Istream& is)
:
    timeIntegratorCoeffs(1),
    order_(readLabel(is)),
    currOrder_(1),
    times_({1.0, 0.0}),
    uniform_(true),
    startIndex_(-1),

    eqn_(times_),
    integrator_(scalarIntegrator::New("Simpson38", eqn_, dictionary()))
{
    if (order_ < 1)
    {
        WarningInFunction
            << "Adams-Bashforth only supports a minimum of 1st order accuracy."
            << endl;
        order_ = 1;
    }
    // else if (order_ > 5)
    // {
    //     WarningInFunction
    //         << "Adams-Bashforth only supports a maximum of 5th order accuracy."
    //         << endl;
    //     order_ = 5;
    // }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeIntegrators::AdamsBashforth::~AdamsBashforth()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeIntegrators::AdamsBashforth::set
(
    List<List<scalar>>& as,
    List<List<scalar>>& bs,
    const label index
) const
{
    if (startIndex_ < 0)
    {
        startIndex_ = index;
    }
    label currOrder = max(min(index - startIndex_, order_), 1);

    if (currOrder < 2)
    {
        as = {{1.0}};
        bs = {{1.0}};
    }
    else if (uniform_ && currOrder < 6)
    {
        if (currOrder == 2)
        {
            as = {{1.0, 0.0}};
            bs = {{1.5, -0.5}};
        }
        else if (currOrder == 3)
        {
            scalar s = 1.0/12.0;
            as = {{1.0, 0.0, 0.0}};
            bs = {{23.0*s, -16.0*s, 5.0*s}};
        }
        else if (currOrder == 4)
        {
            scalar s = 1.0/24.0;
            as = {{1.0, 0.0, 0.0, 0.0}};
            bs = {{55.0*s, -59.0*s, 37.0*s, -9.0*s}};
        }
        else if (currOrder == 5)
        {
            scalar s = 1.0/720.0;
            as = {{1.0, 0.0, 0.0, 0.0, 0.0}};
            bs = {{1901.0*s, -2774.0*s, 2616.0*s, -1274.0*s, 251.0*s}};
        }
    }
    else
    {
        as.setSize(1);
        bs.setSize(1);

        as[0] = List<scalar>(currOrder_, 0.0);
        as[0][0] = 1.0;

        const scalar t = times_[0];
        const scalar t0 = times_[1];
        const scalar dt = t - t0;
        List<scalar> betas(currOrder_);
        forAll(betas, i)
        {
            eqn_.setIndex(i+1);
            betas[i] = integrator_->integrate(t0, t, -1)/dt;
        }

        bs[0].transfer(betas);
    }

}


bool Foam::timeIntegrators::AdamsBashforth::update
(
    const label index,
    const scalar t
)
{
    if(startIndex_ < 0)
    {
        startIndex_ = index;
    }
    currOrder_ = min(index - startIndex_, order_);

    const scalar dt = t - times_[0];

    times_.setSize(currOrder_+1);
    uniform_ = true;
    for (label i = times_.size()-1; i >= 1; i--)
    {
        if (mag(dt - (times_[i-1] - times_[i])) > small)
        {
            uniform_ = false;
        }
        times_[i] = times_[i-1];
    }
    times_[0] = t;

    if  (currOrder_ < order_ || !uniform_)
    {
        return true;
    }
    return false;
}


Foam::List<Foam::label> Foam::timeIntegrators::AdamsBashforth::oldIs
(
    const List<List<scalar>>& as
) const
{
    return List<label>(currOrder_, -1);
}


Foam::List<Foam::label> Foam::timeIntegrators::AdamsBashforth::deltaIs
(
    const List<List<scalar>>& bs
) const
{
    // Only save the first entry
    return identityMap(currOrder_);
}


void Foam::timeIntegrators::AdamsBashforth::setTimeFactors
(
    const List<List<scalar>>& as,
    const List<List<scalar>>& bs,
    List<scalar>& f0,
    List<scalar>& f
) const
{
    f.setSize(nSteps_, 1.0);
    f0.setSize(nSteps_, 0.0);

    // f0[0] = 0.0;
    // f[0] = sum(bs[0]);
    // for (label stepi = 1; stepi < as.size(); stepi++)
    // {
    //     scalarList ts(stepi+1, 0.0);
    //     scalarList dts(stepi, 0.0);
    //     forAll(dts, i)
    //     {
    //         dts[i] = sum(bs[i]);
    //     }
    //     ts[1] = dts[0];
    //
    //     for (label i = 1; i < stepi; i++)
    //     {
    //         for (label j = 0; j < as[i].size(); j++)
    //         {
    //             ts[i+1] += as[i][j]*ts[j];
    //         }
    //         ts[i+1] += dts[i];
    //     }
    //     f0[stepi-1] = ts.last() - dts.last();
    //     f[stepi-1] = f0[stepi-1] + sum(bs[stepi-1]);
    // }
}


Foam::List<Foam::label> Foam::timeIntegrators::AdamsBashforth::deltaSaveMap() const
{
    labelList map(currOrder_, 0);
    for (label i = 0; i < currOrder_-1; i++)
    {
        map[i] = i+1;
    }
    return map;
}

// ************************************************************************* //
