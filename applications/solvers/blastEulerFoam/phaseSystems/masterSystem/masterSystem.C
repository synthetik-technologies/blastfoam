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
#include "packingLimitModel.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(masterSystem, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::masterSystem::masterSystem
(
    const word& type,
    const phaseSystem& fluid,
    const dictionary& dict
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
    dict_(dict),
    group_(dict_.lookupOrDefault<word>("name", type + "Total")),
    writeTotal_(dict_.lookupOrDefault("writeTotal", false)),
    fluid_(fluid),
    phases_(0),
    alphaPtr_(nullptr),
    rhoPtr_(nullptr),
    UPtr_(nullptr),
    phiPtr_(nullptr),
    alphaPhiPtr_(nullptr),
    readResidualAlpha_(dict.found("residualAlpha")),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        dict.lookupOrDefault("residualAlpha", small)
    ),
    readResidualRho_(dict.found("residualRho")),
    residualRho_
    (
        "residualRho",
        dimDensity,
        dict.lookupOrDefault("residualRho", small)
    ),
    alphaMax_
    (
        IOobject
        (
            IOobject::groupName("alphaMax", group_),
            fluid.mesh().time().name(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar(dimless, 1.0),
        extrapolatedCalculatedFvPatchScalarField::typeName
    )
{
    if (writeTotal_)
    {
        this->writeOpt() = IOobject::AUTO_WRITE;
    }
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
                    fluid_.mesh().time().name(),
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
                    fluid_.mesh().time().name(),
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
                    fluid_.mesh().time().name(),
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


Foam::scalar Foam::masterSystem::minAlphaMax() const
{
    return packingLimitModel_->minAlphaMax();
}


void Foam::masterSystem::initialize()
{
    if (!readResidualAlpha_)
    {
        residualAlpha_ = phases_[0].residualAlpha();
        for (label phasei = 1; phasei < phases_.size(); phasei++)
        {
            residualAlpha_ =
                max(residualAlpha_, phases_[phasei].residualAlpha());
        }
    }
    if (!readResidualRho_)
    {
        residualRho_ = phases_[0].residualRho();
        for (label phasei = 1; phasei < phases_.size(); phasei++)
        {
            residualRho_ =
                max(residualRho_, phases_[phasei].residualRho());
        }
    }

    packingLimitModel_ = packingLimitModel::New(dict_, *this);
    packingLimitModel_->updateAlphaMax(alphaMax_);
    alphaMax_.correctBoundaryConditions();

    // Print granular quantities only if more than 1 phase is present
    if (phases_.size() > 1 && !alphaPtr_.valid())
    {
        alphaMax_.writeOpt() = this->writeOpt();

        alphaPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("alpha", group_),
                    fluid_.mesh().time().name(),
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
                    fluid_.mesh().time().name(),
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
    if (phases_.size() > 1)
    {
        alphaMax_.writeOpt() = this->writeOpt();
    }
}


bool Foam::masterSystem::contains(const phaseModel& phase) const
{
    return whichPhase(phase) >= 0;
}


bool Foam::masterSystem::contains(const word& phaseName) const
{
    return whichPhase(phaseName) >= 0;
}


Foam::label Foam::masterSystem::whichPhase(const phaseModel& phase) const
{
    forAll(phases_, phasei)
    {
        if (&phases_[phasei] == &phase)
        {
            return phasei;
        }
    }
    return -1;
}


Foam::label Foam::masterSystem::whichPhase(const word& phaseName) const
{
    forAll(phases_, phasei)
    {
        if (phases_[phasei].name() == phaseName)
        {
            return phasei;
        }
    }
    return -1;
}

void Foam::masterSystem::update()
{
    correctAlpha();

    //- Update packing limit
    if (phases_.size() > 1)
    {
        packingLimitModel_->updateAlphaMax(alphaMax_);
        alphaMax_.correctBoundaryConditions();
    }

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
    for (label phasei = 1; phasei < phases_.size(); phasei++)
    {
        alphaPtr_() += phases_[phasei];
    }
    alphaPtr_().correctBoundaryConditions();
}
// ************************************************************************* //
