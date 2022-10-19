/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "masterSystem.H"
#include "masterSystemList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(masterSystem, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::masterSystem::masterSystem
(
    const word& type,
    const word& group,
    const phaseSystem& fluid
)
:
    regIOobject
    (
        IOobject
        (
            type,
            fluid.mesh().time().constant(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        )
    ),
    group_(group),
    fluid_(fluid),
    phases_(0),
    alphaPtr_(nullptr),
    rhoPtr_(nullptr),
    UPtr_(nullptr),
    phiPtr_(nullptr),
    alphaPhiPtr_(nullptr),
    residualAlpha_("residualAlpha", dimless, small)
{
    masterSystemList::New(fluid.mesh()).addSystem(*this);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::masterSystem::~masterSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volScalarField& Foam::masterSystem::alpha() const
{
    return alphaPtr_.valid() ? alphaPtr_() : phases_[0];
}


const Foam::volScalarField& Foam::masterSystem::rho() const
{
    if (phases_.size() == 1)
    {
        return phases_[0].rho();
    }
    if (!rhoPtr_.valid())
    {
        rhoPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("rho", group_),
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh()
                ),
                phases_[0].rho()
            )
        );
        for (label phasei = 1; phasei < phases_.size(); phasei++)
        {
            rhoPtr_() += phases_[phasei].alphaRho();
        }
        rhoPtr_() /= max(alpha(), residualAlpha_);
    }
    return rhoPtr_();
}


const Foam::volVectorField& Foam::masterSystem::U() const
{
    return UPtr_.valid() ? UPtr_() : phases_[0].U();
}



const Foam::surfaceScalarField& Foam::masterSystem::phi() const
{
    if (phases_.size() == 1)
    {
        return phases_[0].phi();
    }
    if (!phiPtr_.valid())
    {
        phiPtr_.set
        (
            new surfaceScalarField
            (
                IOobject
                (
                    IOobject::groupName("phi", group_),
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh()
                ),
                phases_[0].alphaPhi()
            )
        );
        for (label phasei = 1; phasei < phases_.size(); phasei++)
        {
            phiPtr_() += phases_[phasei].alphaPhi();
        }
        phiPtr_() /= max(fvc::interpolate(alpha()), residualAlpha_);
    }
    return phiPtr_();
}


Foam::tmp<Foam::surfaceScalarField> Foam::masterSystem::alphaPhi() const
{
    if (phases_.size() == 1)
    {
        return phases_[0].alphaPhi();
    }
    if (!alphaPhiPtr_.valid())
    {
        alphaPhiPtr_.set
        (
            new surfaceScalarField
            (
                IOobject
                (
                    IOobject::groupName("alphaPhi", group_),
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh()
                ),
                phases_[0].alphaPhi()
            )
        );
        for (label phasei = 1; phasei < phases_.size(); phasei++)
        {
            alphaPhiPtr_() += phases_[phasei].alphaPhi();
        }
    }
    return alphaPhiPtr_();
}


const Foam::labelList& Foam::masterSystem::phaseIndexes() const
{
    return phaseIndexes_;
}


void Foam::masterSystem::addPhase
(
    phaseModel& phase
)
{
    const label phasei = phases_.size();
    phases_.resize(phasei + 1);
    phases_.set(phasei, &phase);
    phaseIndexes_.append(phase.index());

    // Print granular quantities only if more than 1 phase is present
    if (phases_.size() > 1 && !alphaPtr_.valid())
    {
        alphaPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("alpha", group_),
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh(),
                    IOobject::NO_READ,
                    this->writeOpt()
                ),
                fluid_.mesh(),
                dimensionedScalar("0", dimless, 0.0)
            )
        );
        UPtr_.set
        (
            new volVectorField
            (
                IOobject
                (
                    IOobject::groupName("U", group_),
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh(),
                    IOobject::NO_READ,
                    this->writeOpt()
                ),
                fluid_.mesh(),
                dimensionedVector("0", dimVelocity, Zero)
            )
        );
    }
}


bool Foam::masterSystem::contains(const phaseModel& phase) const
{
    forAll(phases_, phasei)
    {
        if (&phases_[phasei] == &phase)
        {
            return true;
        }
    }
    return false;
}


bool Foam::masterSystem::contains(const word& phaseName) const
{
    forAll(phases_, phasei)
    {
        if (phases_[phasei].name() == phaseName)
        {
            return true;
        }
    }
    return false;
}


Foam::tmp<Foam::volScalarField> Foam::masterSystem::alphaMax() const
{
    scalar minAlphaMax = 1.0;
    forAll(phases_, phasei)
    {
        minAlphaMax = min(minAlphaMax, phases_[phasei].alphaMax());
    }
    return volScalarField::New
    (
        IOobject::groupName("alphaMax", group_),
        fluid_.mesh(),
        minAlphaMax
    );
}

void Foam::masterSystem::update()
{
    correctAlpha();
    if (UPtr_.valid())
    {
        volVectorField& U = UPtr_();
        U = phases_[0]*phases_[0].U();

        for (label phasei = 1; phasei < phases_.size(); phasei++)
        {
            U += phases_[phasei]*phases_[phasei].U();
        }
        U /= max(alpha(), residualAlpha_);
    }

    if (rhoPtr_.valid())
    {
        volScalarField& rho = rhoPtr_();
        rho = phases_[0].alphaRho();
        for (label phasei = 1; phasei < phases_.size(); phasei++)
        {
            rho += phases_[phasei].alphaRho();
        }
        rho /= max(alpha(), residualAlpha_);
    }

    if (alphaPhiPtr_.valid())
    {
        surfaceScalarField& alphaPhi = alphaPhiPtr_();
        alphaPhi = phases_[0].alphaPhi();
        for (label phasei = 1; phasei < phases_.size(); phasei++)
        {
            alphaPhi += phases_[phasei].alphaPhi();
        }
    }

    if (phiPtr_.valid())
    {
        phiPtr_() = alphaPhi()/max(fvc::interpolate(alpha()), residualAlpha_);
    }

}

void Foam::masterSystem::correctAlpha()
{
    if (!alphaPtr_.valid())
    {
        return;
    }

    alphaPtr_() = phases_[0];
    for (label phasei = 0; phasei < phases_.size(); phasei++)
    {
        alphaPtr_() += phases_[phasei];
    }
    alphaPtr_().correctBoundaryConditions();
}
// ************************************************************************* //
