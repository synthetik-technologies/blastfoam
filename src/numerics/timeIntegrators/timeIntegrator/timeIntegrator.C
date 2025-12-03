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

#include "timeIntegrator.H"
#include "timeIntegrationSystemBase.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(timeIntegrator, 0);
}


// * * * * * * * * * * * * * * Protected Functions * * * * * * * * * * * * * //


void Foam::timeIntegrator::updateCoeffs()
{
    if (coeffs_->update(curTimeIndex_, time().value()))
    {
        initialize();
    }
    forAll(systems_, i)
    {
        systems_[i].save();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeIntegrator::timeIntegrator
(
    const objectRegistry& obr,
    const dictionary& dict,
    const bool embedded
)
:
    regIOobject
    (
        IOobject
        (
            typeName,
            obr.time().name(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    obr_(obr),
    coeffs_(timeIntegratorCoeffs::New(dict, embedded)),
    systems_(0),
    nSteps_(0),
    stepi_(0),
    as_(0),
    bs_(0),
    f_(0),
    f0_(0),
    oldIs_(0),
    nOld_(0),
    deltaIs_(0),
    nDelta_(0),
    curTimeIndex_(-1),
    restart_(false)
{
    if (!embedded)
    {
        initialize();
    }
}

Foam::timeIntegrator::timeIntegrator
(
    const objectRegistry& obr,
    const dictionary& dict
)
:
    timeIntegrator(obr, dict, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeIntegrator::~timeIntegrator()
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

void Foam::timeIntegrator::initialize()
{
    coeffs_->set(as_, bs_, time().timeIndex());
    nSteps_ = coeffs_->nSteps();

    as_.setSize(nSteps_);
    bs_.setSize(nSteps_);

    oldIs_ = coeffs_->oldIs(as_);
    nOld_ = 0;
    forAll(oldIs_, i)
    {
        nOld_ += oldIs_[i] >= 0;
    }

    deltaIs_ = coeffs_->deltaIs(bs_);
    nDelta_ = 0;
    forAll(deltaIs_, i)
    {
        nDelta_ += deltaIs_[i] >= 0;
    }

    coeffs_->setTimeFactors(as_, bs_, f0_, f_);
}

void Foam::timeIntegrator::addSystem(timeIntegrationSystemBase& system)
{
    DebugInfo<< "Adding timeIntegrationSystem " << system.name() << endl;
    system.set(*this);
    label oldSize = systems_.size();
    systems_.resize(oldSize + 1);
    systems_.set(oldSize, &system);
}


void Foam::timeIntegrator::update()
{
    forAll(systems_, i)
    {
        systems_[i].update();
    }
}


void Foam::timeIntegrator::preUpdate()
{
    if (obr_.time().subCycling())
    {
        curTimeIndex_ = obr_.time().timeIndex();
        restart_ = false;
        updateCoeffs();
    }
    else if
    (
        obr_.time().timeIndex() == curTimeIndex_
     && !obr_.time().subCycling()
    )
    {
        reset();
        restart_ = true;
        Info<< "Restarting time step" << endl;
    }
    else
    {
        curTimeIndex_ = obr_.time().timeIndex();
        restart_ = false;

        updateCoeffs();
    }

    forAll(systems_, i)
    {
        systems_[i].preUpdate();
    }
}


void Foam::timeIntegrator::integrate()
{
    integrate(true, true, true, true, true);
}


void Foam::timeIntegrator::integrate
(
    const bool doExplicit,
    const bool doStore,
    const bool doImplicit,
    const bool doPost,
    const bool doClear
)
{
    preUpdate();

    // Update and store original fields
    for (stepi_ = 0; stepi_ < coeffs_->nSteps(); stepi_++)
    {
        Info<< coeffs_->type() << ": step " << stepi_ << endl;
        this->update();
        forAll(systems_, i)
        {
            Info<< "Solving " << systems_[i].name() << ":" << endl;
            systems_[i].solve();
            Info<< endl;
        }
    }
    stepi_ = coeffs_->nSteps()-1;

    stepi_ = -1;

    if (doExplicit)
    {
        solveExplicit();
    }

    if (doStore)
    {
        storeExplicit();
    }

    if (doImplicit)
    {
        solveImplicit();
    }

    if (doPost)
    {
        postUpdate();
    }

    if (doClear)
    {
        clear();
    }
}


void Foam::timeIntegrator::solveExplicit()
{
    forAll(systems_, i)
    {
        systems_[i].solveExplicit();
    }
}


void Foam::timeIntegrator::storeExplicit()
{
    forAll(systems_, i)
    {
        systems_[i].storeExplicit();
    }
}


void Foam::timeIntegrator::solveImplicit()
{
    forAll(systems_, i)
    {
        systems_[i].solveImplicit();
    }
}


void Foam::timeIntegrator::postUpdate()
{
    forAll(systems_, i)
    {
        systems_[i].postUpdate();
    }
}


Foam::tmp<Foam::scalarField> Foam::timeIntegrator::V0() const
{
    return tmp<scalarField>();
}


Foam::tmp<Foam::scalarField> Foam::timeIntegrator::V() const
{
    return tmp<scalarField>();
}


Foam::scalar Foam::timeIntegrator::totalV0() const
{
    return 1.0;
}


Foam::scalar Foam::timeIntegrator::totalV() const
{
    return 1.0;
}


void Foam::timeIntegrator::clear()
{
    if (!coeffs_->save())
    {
        forAll(systems_, i)
        {
            systems_[i].clear();
        }
    }
}


void Foam::timeIntegrator::reset()
{
    forAll(systems_, i)
    {
        systems_[i].reset();
    }
}

// ************************************************************************* //
