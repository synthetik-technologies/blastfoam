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

#include "kineticTheorySystem.H"
#include "kineticTheoryModel.H"
#include "packingLimitModel.H"
#include "radialModel.H"
#include "viscosityModel.H"
#include "frictionalStressModel.H"
#include "granularPressureModel.H"
#include "conductivityModel.H"
#include "phaseSystem.H"
#include "mathematicalConstants.H"
#include "SortableList.H"
#include "zeroGradientFvPatchFields.H"
#include "noneViscosity.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(kineticTheorySystem, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheorySystem::kineticTheorySystem
(
    const phaseSystem& fluid
)
:
    masterSystem
    (
        typeName,
        fluid,
        fluid.subDict("kineticTheory")
    ),
    ThetapPtr_(nullptr),
    kineticTheoryModels_(0),
    Thetas_(0),
    radialModel_
    (
        kineticTheoryModels::radialModel::New
        (
            dict_,
            *this
        )
    ),
    viscosityModel_
    (
        kineticTheoryModels::viscosityModel::New
        (
            dict_,
            *this
        )
    ),
    granularPressureModel_
    (
        kineticTheoryModels::granularPressureModel::New
        (
            dict_,
            *this
        )
    ),
    conductivityModel_
    (
        kineticTheoryModels::conductivityModel::New
        (
            dict_,
            *this
        )
    ),
    frictionalStressModel_
    (
        kineticTheoryModels::frictionalStressModel::New(dict_, *this)
    ),
    eTable_(dict_.lookupOrDefault("e", phasePair::scalarTable())),
    CfTable_(dict_.lookupOrDefault("Cf", phasePair::scalarTable())),
    alphaMinFriction_
    (
        IOobject
        (
            IOobject::groupName("alphaMinFriction", group_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("one", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    minTheta_
    (
        dimensionedScalar::lookupOrDefault
        (
            "minTheta",
            dict_,
            sqr(dimVelocity),
            1e-8
        )
    ),
    includeViscosity_
    (
        !isA<kineticTheoryModels::noneViscosity>(viscosityModel_())
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheorySystem::~kineticTheorySystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::kineticTheorySystem::es(const phasePairKey& pair) const
{
    label i = phaseIndexes_[fluid_.phases()[pair[0]].index()];
    if (pair[0] == pair[1])
    {
        return kineticTheoryModels_[i].es();
    }
    if (eTable_.found(pair))
    {
        return eTable_[pair];
    }

    label j = phaseIndexes_[fluid_.phases()[pair[1]].index()];
    return
        (
            kineticTheoryModels_[i].es()
          + kineticTheoryModels_[j].es()
        )/2.0;
}

Foam::scalar Foam::kineticTheorySystem::Cf(const phasePairKey& pair) const
{
    if (CfTable_.found(pair))
    {
        return CfTable_[pair];
    }
    return 0.0;
}


Foam::tmp<Foam::volScalarField> Foam::kineticTheorySystem::gs0
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const bool calc
) const
{
    if (&phase1  == &phase2 && !calc)
    {
        return kineticTheoryModels_[whichPhase(phase1)].gs0();
    }
    return radialModel_->gs0(phase1, phase2);
}


Foam::scalar Foam::kineticTheorySystem::cellgs0
(
    const label celli,
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    if (&phase1  == &phase2)
    {
        return kineticTheoryModels_[whichPhase(phase1)].gs0()[celli];
    }
    return radialModel_->cellgs0(celli, phase1, phase2);
}


Foam::tmp<Foam::volScalarField> Foam::kineticTheorySystem::gs0Prime
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const bool calc
) const
{
    if (&phase1  == &phase2 && !calc)
    {
        return kineticTheoryModels_[whichPhase(phase1)].gs0Prime();
    }
    return radialModel_->gs0prime(phase1, phase2);
}


Foam::scalar Foam::kineticTheorySystem::cellgs0Prime
(
    const label celli,
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    if (&phase1  == &phase2)
    {
        return kineticTheoryModels_[whichPhase(phase1)].gs0Prime()[celli];
    }
    return radialModel_->cellgs0prime(celli, phase1, phase2);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheorySystem::nu
(
    const phaseModel& phase,
    const volScalarField& Theta
) const
{
    phasePairKey key(phase.name(), phase.name(), false);
    return viscosityModel_->nu
    (
        phase,
        Theta,
        kineticTheoryModels_[whichPhase(phase)].gs0(),
        phase.rho(),
        phase.d(),
        dimensionedScalar("e", dimless, es(key))
    );
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheorySystem::Ps(const phaseModel& phase) const
{
    tmp<volScalarField> tmpPs
    (
        new volScalarField
        (
            IOobject::groupName("Ps", phase.name()),
            phase*phase.rho()*phase.Theta()
        )
    );
    volScalarField& ps = tmpPs.ref();

    forAll(phases_, phasej)
    {
        const phaseModel& phase2 = phases_[phasej];
        phasePairKey key(phase.name(), phase2.name(), false);

        ps += granularPressureModel_->granularPressure
        (
            phase,
            phase2,
            phase.Theta(),
            phase2.Theta(),
            gs0(phase, phase2),
            es(key)
        );
    }
    return tmpPs;
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheorySystem::Ps
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    phasePairKey key(phase1.name(), phase2.name(), false);
    return granularPressureModel_->granularPressure
    (
        phase1,
        phase2,
        phase1.Theta(),
        phase2.Theta(),
        gs0(phase1, phase2),
        es(key)
    );
}



Foam::tmp<Foam::volScalarField>
Foam::kineticTheorySystem::dPsdAlpha(const phaseModel& phase) const
{
    tmp<volScalarField> tmpdPsdAlpha
    (
        new volScalarField
        (
            IOobject::groupName("dPsdAlpha", phase.name()),
            phase.rho()*phase.Theta()
        )
    );
    volScalarField& dPsdAlpha = tmpdPsdAlpha.ref();

    forAll(phaseIndexes_, phasej)
    {
        const phaseModel& phase2 = phases_[phasej];
        phasePairKey key(phase.name(), phase2.name(), false);

        dPsdAlpha += granularPressureModel_->granularPressureByAlpha
        (
            phase,
            phase2,
            phase.Theta(),
            phase2.Theta(),
            gs0(phase, phase2),
            gs0Prime(phase, phase2),
            es(key)
        );
    }
    return tmpdPsdAlpha;
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheorySystem::dPsdTheta(const phaseModel& phase) const
{
    tmp<volScalarField> tmpdPsdTheta
    (
        new volScalarField
        (
            IOobject::groupName("dPsdTheta", phase.name()),
            phase*phase.rho()
        )
    );
    volScalarField& dPsdTheta = tmpdPsdTheta.ref();

    forAll(phases_, phasej)
    {
        const phaseModel& phase2 = phases_[phasej];
        phasePairKey key(phase.name(), phase2.name(), false);

        dPsdTheta += granularPressureModel_->granularPressureByTheta
        (
            phase,
            phase2,
            phase.Theta(),
            phase2.Theta(),
            gs0(phase, phase2),
            es(key)
        );
    }
    return tmpdPsdTheta;
}

Foam::tmp<Foam::volScalarField>
Foam::kineticTheorySystem::kappa
(
    const phaseModel& phase,
    const volScalarField& Theta
) const
{
    phasePairKey key(phase.name(), phase.name(), false);
    return conductivityModel_->kappa
    (
        phase,
        Theta,
        kineticTheoryModels_[phaseIndexes_[phase.index()]].gs0(),
        phase.rho(),
        phase.d(),
        dimensionedScalar("e", dimless, es(key))
    );
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheorySystem::lambda
(
    const phaseModel& phase
) const
{
    // Bulk viscosity as a function of all phases (Eq. 10, pg. 3780)
    tmp<volScalarField> tmpLambda
    (
        new volScalarField
        (
            IOobject
            (
                "lambda",
                fluid_.mesh().time().timeName(),
                fluid_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fluid_.mesh(),
            dimensionedScalar(dimensionSet(0, 2, -1, 0, 0), 0)
        )
    );
    volScalarField& l = tmpLambda.ref();

    using Foam::constant::mathematical::pi;
    tmp<volScalarField> tm1(phase.rho()*pi*pow3(phase.d())/6.0);
    const volScalarField& m1 = tm1();

    tmp<volScalarField> tTheta1 = phase.Theta();
    const volScalarField& Theta1 = tTheta1();

    forAll(phases_, phasej)
    {
        const phaseModel& phase2 = phases_[phasej];

        tmp<volScalarField> tm2(phase2.rho()*pi*pow3(phase2.d())/6.0);
        const volScalarField&  m2 = tm2();

        tmp<volScalarField> tTheta2 = phase2.Theta();
        const volScalarField& Theta2 = tTheta2();

        tmp<volScalarField> tPsij(Ps(phase, phase2));
        const volScalarField& Psij = tPsij();

        l +=
            Psij/phase.rho()*(phase.d() + phase2.d())/6.0
           *sqrt
            (
                2.0*sqr(m1*Theta1 + m2*Theta2)
               /max
                (
                    pi*Theta1*Theta2*(sqr(m1)*Theta1 + sqr(m2)*Theta2),
                    dimensionedScalar(dimensionSet(2, 6, -6, 0, 0), small)
                )
            );
    }
    return tmpLambda;
}


Foam::tmp<Foam::volScalarField> Foam::kineticTheorySystem::frictionalPressure
(
    const phaseModel& phase
) const
{
    if (this->polydisperse())
    {
        return
            phase/max(this->alpha(), this->residualAlpha())
           *frictionalStressModel_->frictionalPressure
            (
                phase,
                this->alpha(),
                alphaMax_
            );
    }
    else
    {
        return frictionalStressModel_->frictionalPressure
        (
            phase,
            this->alpha(),
            alphaMax_
        );
    }
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheorySystem::frictionalPressurePrime
(
    const phaseModel& phase
) const
{
    if (this->polydisperse())
    {
        return
        (
            (this->alpha() - phase)
           *frictionalStressModel_->frictionalPressure
            (
                phase,
                this->alpha(),
                alphaMax_
            )
          + phase*this->alpha()
           *frictionalStressModel_->frictionalPressurePrime
            (
                phase,
                this->alpha(),
                alphaMax_
            )
        )/sqr(max(this->alpha(), this->residualAlpha()));
    }
    else
    {
        return frictionalStressModel_->frictionalPressurePrime
        (
            phase,
            this->alpha(),
            alphaMax_
        );
    }
}


Foam::tmp<Foam::volScalarField> Foam::kineticTheorySystem::muFrictional
(
    const phaseModel& phase,
    const volScalarField& Pfr
) const
{
    return frictionalStressModel_->mu
    (
        phase,
        this->alpha(),
        alphaMax_,
        Pfr
    );
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheorySystem::productionSource
(
    const kineticTheoryModel& kt,
    const phaseModel& phase2
) const
{
    const phaseModel& phase1 = kt.phase();

    // Production of granular energy (Houim and Oran 2016, Eq. 3.49, Eq. B 66)
    return tmp<volScalarField>
    (
        new volScalarField
        (
            81.0*phase1*sqr(phase2.mu())*magSqr(phase1.U() - phase2.U())
           /(
                kt.gs0()*pow3(phase1.d())
               *phase1.rho()
               *sqrt(Foam::constant::mathematical::pi)
            )
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheorySystem::dissipationSource
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const dimensionedScalar& deltaT
) const
{
    const scalar pi(Foam::constant::mathematical::pi);
    if (Thetas_.size() == 1)
    {
        const kineticTheoryModel& kt = kineticTheoryModels_[0];
        if (kt.es() == 1)
        {
            return volScalarField::New
            (
                "dissipationSource." + phase1.group(),
                phase1.mesh(),
                dimensionedScalar(dimDensity*sqr(dimVelocity), 0.0)
            );
        }
        tmp<volScalarField> gammaCoeff
        (
            volScalarField::New
            (
                "gammaCoeff",
                12.0
               *(1.0 - sqr(kt.es()))
               *kt.gs0()
               *phase1*phase1.alphaRho()
               /(phase1.d()*sqrt(pi))
            )
        );
        tmp<volScalarField> ThetaStar
        (
            kt.Theta()
           *sqr
            (
                3.0*phase1.alphaRho()
               /(
                    max
                    (
                        3.0*phase1.alphaRho()
                      + deltaT*gammaCoeff*sqrt(kt.Theta()),
                        3.0*phase1.residualAlphaRho()
                    )
                )
            )
        );
        return volScalarField::New
        (
            "dissipationSource." + phase1.group(),
            1.5*phase1.alphaRho()*(ThetaStar - kt.Theta())
        );
    }

    phasePairKey key(phase1.name(), phase2.name(), false);
    const scalar e = this->es(key);
    if (e == 1)
    {
        return volScalarField::New
        (
            "dissipationSource." + phase1.group() + "." + phase2.group(),
            phase1.mesh(),
            dimensionedScalar(dimDensity*sqr(dimVelocity), 0.0)
        );
    }

    // Dissipation of granular energy (Huilin and Gidaspow 2003, Eq. 25)
    volScalarField Theta1(phase1.Theta());
    Theta1.max(1e-10);
    volScalarField Theta2(phase2.Theta());
    Theta2.max(1e-10);

    tmp<volScalarField> m1(pi/6.0*pow3(phase1.d())*phase1.rho());
    tmp<volScalarField> m2(pi/6.0*pow3(phase2.d())*phase2.rho());
    tmp<volScalarField> m0(m1() + m2());
    tmp<volScalarField> m1Thetam2Theta(sqr(m1())*Theta1 + sqr(m2())*Theta2);

    return volScalarField::New
    (
        "dissipationSource." + phase1.group() + "." + phase2.group(),
       - (
            (
                3.0/phase1.d()
               *sqrt
                (
                    2.0*sqr(m0())*phase1.Theta()*phase2.Theta()
                   /(pi*m1Thetam2Theta())
                )
              - (3.0*m0()*(m1()*phase1.Theta() + m2()*phase2.Theta()))
               /(4.0*m1Thetam2Theta())
               *fvc::div(phase1.phi())
            )
           *(1.0 - e)
           *this->Ps(phase1, phase2)
        )*deltaT
    );
}


void Foam::kineticTheorySystem::addPhase
(
    phaseModel& phase
)
{
    const label phasei = phases_.size();
    masterSystem::addPhase(phase);
    kineticTheoryModels_.resize(phasei + 1);
    Thetas_.resize(phasei + 1);

    kineticTheoryModel& kt = dynamicCast<kineticTheoryModel>(phase);
    kineticTheoryModels_.set(phasei, &kt);
    Thetas_.set(phasei, &kt.Theta());
}

void Foam::kineticTheorySystem::initialize()
{
    masterSystem::initialize();

    // Print granular quantities only if more than 1 phase is present
    if (phases_.size() > 1 && !ThetapPtr_.valid())
    {
        ThetapPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("Theta", group_),
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh(),
                    IOobject::NO_READ,
                    this->writeOpt()
                ),
                fluid_.mesh(),
                dimensionedScalar("0", sqr(dimVelocity), 0.0)
            )
        );
    }

    forAll(phases_, phasei)
    {
        word phaseName = phases_[phasei].name();
        forAll(phases_, phasej)
        {
            word otherPhaseName = phases_[phasej].name();
            phasePairKey key
            (
                phaseName,
                otherPhaseName,
                false
            );
            pairs_.append(key);
            word name(key.second());
            name[0] = toupper(name[0]);
            name = key.first() + "And" + name;

            if (phaseName == otherPhaseName)
            {
                name = phaseName;
            }
        }
    }

    update();
}

void Foam::kineticTheorySystem::update()
{
    masterSystem::update();
    if (ThetapPtr_.valid())
    {
        volScalarField& Thetap = ThetapPtr_();
        Thetap = phases_[0]*Thetas_[0];
        for (label phasei = 1; phasei < phases_.size(); phasei++)
        {
            Thetap +=  phases_[phasei]*Thetas_[phasei];
        }
        Thetap /= max(alpha(), residualAlpha_);
    }

    frictionalStressModel_->update();
    alphaMinFriction_ =
        frictionalStressModel_->alphaMinFriction(this->alpha(), alphaMax_);
}


void Foam::kineticTheorySystem::solve()
{
    frictionalStressModel_->solve();
}


void Foam::kineticTheorySystem::postUpdate()
{
    frictionalStressModel_->postUpdate();
}


void Foam::kineticTheorySystem::clear()
{
    frictionalStressModel_->clear();
}

// ************************************************************************* //
