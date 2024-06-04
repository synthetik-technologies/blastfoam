/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
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

#include "embeddedTimeIntegrator.H"
#include "embeddedTimeIntegrationSystemBase.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(embeddedTimeIntegrator, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::embeddedTimeIntegrator::embeddedTimeIntegrator
(
    const objectRegistry& obr,
    const dictionary& dict
)
:
    timeIntegrator(obr, dict, true),
    dict_(obr.time().controlDict().subOrEmptyDict("adaptiveStepControls")),
    safeScale_(dict_.lookupOrDefault<scalar>("safeScale", 0.9)),
    alphaInc_(dict_.lookupOrDefault<scalar>("alphaIncrease", 0.2)),
    alphaDec_(dict_.lookupOrDefault<scalar>("alphaDecrease", 0.25)),
    minScale_(dict_.lookupOrDefault<scalar>("minScale", 0.2)),
    maxScale_(dict_.lookupOrDefault<scalar>("maxScale", 10)),
    deltaT_(obr.time().deltaTValue()),
    adjust_(true)
{
    initialize();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::embeddedTimeIntegrator::~embeddedTimeIntegrator()
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

void Foam::embeddedTimeIntegrator::initialize()
{
    coeffs_->set(as_, bs_, time().timeIndex());
    nSteps_ = coeffs_->nSteps();

    // Make sure as and bs are the right size
    as_.setSize(nSteps_);

    // Set required old times to save
    oldIs_ = coeffs_->oldIs(as_);
    if (!oldIs_.size())
    {
        oldIs_.append(0);
    }
    else if (oldIs_[0] < 0)
    {
        forAll(oldIs_, i)
        {
            if (i == 0 || oldIs_[i] >= 0)
            {
                oldIs_[i]++;
            }
        }
    }
    nOld_ = 0;
    forAll(oldIs_, i)
    {
        nOld_ += oldIs_[i] >= 0;
    }

    // Set required deltas to save
    deltaIs_ = identityMap(nSteps_);
    nDelta_ = nSteps_;

    // Set error coefficients b_{n-2} - b_{n-1}
    const scalarList& highOrderCoeffs = bs_[bs_.size()-2];
    scalarList& errorCoeffs = bs_[bs_.size()-1];
    forAll(errorCoeffs, i)
    {
        errorCoeffs[i] = highOrderCoeffs[i] - errorCoeffs[i];
    }

    coeffs_->setTimeFactors(as_, bs_, f0_, f_);
}

void Foam::embeddedTimeIntegrator::addSystem
(
    timeIntegrationSystemBase& system
)
{
    DebugInfo<< "Adding timeIntegrationSystem " << system.name() << endl;
    if (!isA<embeddedTimeIntegrationSystemBase>(system))
    {
        FatalErrorInFunction
            << "Adding a system without error calculation to an "
            << "embedded time integrator" << endl
            << abort(FatalError);
    }
    timeIntegrator::addSystem(system);
}


void Foam::embeddedTimeIntegrator::integrate()
{
    const Time& runTime = obr_.time();

//     if (runTime.timeIndex() == curTimeIndex_)
//     {
//         reset();
//         restart_ = true;
//         DebugInfo<< "Restarting time step" << endl;
//     }
//     else
    {
        curTimeIndex_ = obr_.time().timeIndex();
        restart_ = false;

        update();
    }

    const scalar t0 = runTime.value() - runTime.deltaTValue();
    scalar dt = runTime.deltaTValue();

    // Begin ode loop
    scalar err = 0.0;
    do
    {
        // Update and store original fields
        for (stepi_ = 0; stepi_ < coeffs_->nSteps(); stepi_++)
        {
            this->updateAll();
            forAll(systems_, i)
            {
                systems_[i].solve();
            }
        }

        err = 0.0;
        forAll(systems_, i)
        {
            err =
                max
                (
                    err,
                    dynamicCast<embeddedTimeIntegrationSystemBase>
                    (
                        systems_[i]
                    ).error()
                );
        }

        if (err > 1 && adjust_)
        {
            scalar scale =
                max(safeScale_*pow(err, -alphaDec_), minScale_);
            dt *= scale;
            const_cast<Time&>(runTime).setDeltaTNoAdjust(dt);

            reset();

            if (dt < vSmall)
            {
                FatalErrorInFunction
                    << "0 sized time step" << endl
                    << abort(FatalError);
            }
        }
    } while (err > 1);

    if (adjust_)
    {
        const_cast<Time&>(runTime).setTime(t0 + dt, runTime.timeIndex());
        if (err > pow(maxScale_/safeScale_, -1.0/alphaInc_))
        {
            deltaT_ =
                min
                (
                    max
                    (
                        safeScale_*pow(err, -alphaInc_),
                        minScale_
                    ),
                    maxScale_
                )*dt;
        }
        else
        {
            deltaT_ = dt*safeScale_*maxScale_;
        }
    }
    else
    {
        deltaT_ = runTime.deltaTValue();
    }

    this->postUpdateAll();
    stepi_ = -1;
}

// ************************************************************************* //
