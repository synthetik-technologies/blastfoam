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


Foam::kineticTheorySystem& Foam::kineticTheorySystem::New
(
    const phaseSystem& fluid
)
{
    const fvMesh& mesh = fluid.mesh();
    if (!mesh.foundObject<kineticTheorySystem>(kineticTheorySystem::typeName))
    {
        kineticTheorySystem* ktPtr
        (
            new kineticTheorySystem(fluid)
        );

        // Transfer ownership of this object to the objectRegistry
        ktPtr->store(ktPtr);
    }

    return mesh.lookupObjectRef<kineticTheorySystem>(kineticTheorySystem::typeName);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheorySystem::kineticTheorySystem
(
    const phaseSystem& fluid
)
:
    regIOobject
    (
        IOobject
        (
            "kineticTheorySystem",
            fluid.mesh().time().constant(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        )
    ),
    fluid_(fluid),
    dict_(fluid.subDict("kineticTheory")),
    name_
    (
        dict_.lookupOrDefault<word>
        (
            "name",
            "kineticTheoryTotal"
        )
    ),
    writeTotal_(dict_.lookupOrDefault("writeTotal", false)),
    alphapPtr_(nullptr),
    UpPtr_(nullptr),
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
        kineticTheoryModels::frictionalStressModel::New(dict_)
    ),
    eTable_(dict_.lookupOrDefault("e", phasePair::scalarTable())),
    CfTable_(dict_.lookupOrDefault("Cf", phasePair::scalarTable())),
    alphaMax_
    (
        IOobject
        (
            IOobject::groupName("alphaMax", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("one", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    minAlphaMax_(1.0),
    alphaMinFriction_
    (
        IOobject
        (
            IOobject::groupName("alphaMinFriction_", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("one", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    residualAlpha_("residualAlpha", dimless, dict_),
    Pfr_
    (
        IOobject
        (
            IOobject::groupName("Pfr", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("zero", dimPressure, 0.0)
    ),

    PfrPrime_
    (
        IOobject
        (
            IOobject::groupName("PfrPrime", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("zero", dimPressure, 0.0)
    ),

    nuFric_
    (
        IOobject
        (
            IOobject::groupName("nuFric", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar(dimensionSet(0, 2, -1, 0, 0), 0)
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

bool Foam::kineticTheorySystem::read()
{
    residualAlpha_.readIfPresent(dict_);

    radialModel_->read();
    viscosityModel_->read();
    frictionalStressModel_->read();
    granularPressureModel_->read();
    conductivityModel_->read();

    writeTotal_ = dict_.lookupOrDefault("writeTotal", false);

    if (alphapPtr_.valid())
    {
        if (writeTotal_)
        {
            alphapPtr_->writeOpt() = AUTO_WRITE;
            UpPtr_->writeOpt() = AUTO_WRITE;
            ThetapPtr_->writeOpt() = AUTO_WRITE;
            alphaMax_.writeOpt() = AUTO_WRITE;
        }
        else
        {
            alphapPtr_->writeOpt() = NO_WRITE;
            UpPtr_->writeOpt() = NO_WRITE;
            ThetapPtr_->writeOpt() = NO_WRITE;
            alphaMax_.writeOpt() = NO_WRITE;
        }
    }

    return true;
}

bool Foam::kineticTheorySystem::readIfModified()
{
    return true;
}


bool Foam::kineticTheorySystem::polydisperse() const
{
    return (phaseIndexes_.size() > 1);
}


const Foam::volScalarField& Foam::kineticTheorySystem::alphap() const
{
    return alphapPtr_.valid() ? alphapPtr_() : kineticTheoryModels_[0].phase();
}


const Foam::volVectorField& Foam::kineticTheorySystem::Up() const
{
    return UpPtr_.valid() ? UpPtr_() : kineticTheoryModels_[0].phase().U();
}


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
    return CfTable_[pair];
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
        return kineticTheoryModels_[phaseIndexes_[phase1.index()]].gs0();
    }
    return radialModel_->gs0(phase1, phase2);
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
        return kineticTheoryModels_[phaseIndexes_[phase1.index()]].gs0Prime();
    }
    return radialModel_->gs0prime(phase1, phase2);
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
        kineticTheoryModels_[phaseIndexes_[phase.index()]].gs0(),
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

    forAll(ktPhases_, phasej)
    {
        const phaseModel& phase2 = ktPhases_[phasej];
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
        const phaseModel& phase2 = ktPhases_[phasej];
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

    forAll(ktPhases_, phasej)
    {
        const phaseModel& phase2 = ktPhases_[phasej];
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

    scalar pi = Foam::constant::mathematical::pi;
    volScalarField m1(phase.rho()*pi*pow3(phase.d())/6.0);
    const volScalarField& Theta1 = phase.Theta();
    forAll(ktPhases_, phasej)
    {
        const phaseModel& phase2 = ktPhases_[phasej];
        volScalarField m2(phase2.rho()*pi*pow3(phase2.d())/6.0);
        const volScalarField& Theta2 = phase2.Theta();
        volScalarField Psij(Ps(phase, phase2));

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


const Foam::volScalarField& Foam::kineticTheorySystem::frictionalPressure() const
{
    return Pfr_;
}


const Foam::volScalarField&
Foam::kineticTheorySystem::frictionalPressurePrime
(
    const volScalarField& alpha
) const
{
    return PfrPrime_;
}


const Foam::volScalarField& Foam::kineticTheorySystem::nuFrictional() const
{
    return nuFric_;
}


const Foam::labelList& Foam::kineticTheorySystem::phaseIndexes() const
{
    return phaseIndexes_;
}


void Foam::kineticTheorySystem::addPackingPhase(const phaseModel& phase)
{
    label phasei = packingPhases_.size();
    packingPhases_.resize(phasei + 1);
    packingPhases_.set(phasei, &phase);

    minAlphaMax_ = min(minAlphaMax_, phase.alphaMax());

    if (!packingLimitModel_.valid())
    {
        packingLimitModel_ =
        (
            kineticTheoryModels::packingLimitModel::New
            (
                dict_,
                *this
            )
        );
    }

    if (packingPhases_.size() > 1 && !alphapPtr_.valid())
    {
        IOobject::writeOption writeOpt =
            writeTotal_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE;

        alphapPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("alpha", name_),
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh(),
                    IOobject::NO_READ,
                    writeOpt
                ),
                fluid_.mesh(),
                dimensionedScalar("0", dimless, 0.0)
            )
        );
        UpPtr_.set
        (
            new volVectorField
            (
                IOobject
                (
                    IOobject::groupName("U", name_),
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh(),
                    IOobject::NO_READ,
                    writeOpt
                ),
                fluid_.mesh(),
                dimensionedVector("0", dimVelocity, Zero)
            )
        );
        alphaMax_.writeOpt() = writeOpt;
    }
}


void Foam::kineticTheorySystem::addPhase
(
    kineticTheoryModel& kt
)
{
    const phaseModel& phase = kt.phase();
    addPackingPhase(phase);
    word phaseName(phase.name());
    label phasei = kineticTheoryModels_.size();
    kineticTheoryModels_.resize(phasei + 1);
    ktPhases_.resize(phasei + 1);
    Thetas_.resize(phasei + 1);

    ktPhases_.set(phasei, &phase);
    kineticTheoryModels_.set(phasei, &kt);
    Thetas_.set(phasei, &kt.Theta());
    phaseIndexes_.append(phase.index());

    // Print granular quantities only if more than 1 phase is present
    if (ktPhases_.size() > 1 && !ThetapPtr_.valid())
    {
        IOobject::writeOption writeOpt =
            writeTotal_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE;
        ThetapPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("Theta", name_),
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh(),
                    IOobject::NO_READ,
                    writeOpt
                ),
                fluid_.mesh(),
                dimensionedScalar("0", sqr(dimVelocity), 0.0)
            )
        );
    }

    forAll(ktPhases_, phasej)
    {
        word otherPhaseName = ktPhases_[phasej].name();
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


bool Foam::kineticTheorySystem::found(const word& phaseName) const
{
    forAll(ktPhases_, phasei)
    {
        if (ktPhases_[phasei].name() == phaseName)
        {
            return true;
        }
    }
    return false;
}


void Foam::kineticTheorySystem::correct()
{
    if (alphapPtr_.valid())
    {
        volScalarField& alphap = alphapPtr_();
        volVectorField& Up = UpPtr_();

        alphap = Zero;
        Up = Zero;

        forAll(packingPhases_, phasei)
        {
            const phaseModel& phase = packingPhases_[phasei];
            const volScalarField& alpha = phase;

            alphap += alpha;
            Up += alpha*phase.U();
        }
        Up /= max(alphap, residualAlpha_);
    }
    alphaMax_ = packingLimitModel_->alphaMax();
    alphaMax_.correctBoundaryConditions();

    if (!ktPhases_.size())
    {
        return;
    }

    if (ThetapPtr_.valid())
    {
        volScalarField sumAlpha(ktPhases_[0]);
        volScalarField& Thetap = ThetapPtr_();
        Thetap = Thetas_[0];

        for (label phasei = 0; phasei < ktPhases_.size(); phasei++)
        {
            const phaseModel& phase = ktPhases_[phasei];
            const volScalarField& alpha = phase;
            Thetap +=  alpha*Thetas_[phasei];
            sumAlpha += alpha;
        }
        Thetap /= max(sumAlpha, residualAlpha_);
    }

    const volScalarField& alpha = this->alphap();

    alphaMinFriction_ =
        frictionalStressModel_->alphaMinFriction(alpha, alphaMax_);

    Pfr_ = frictionalStressModel_->frictionalPressure
    (
        alpha,
        alphaMax_
    );

    PfrPrime_ = frictionalStressModel_->frictionalPressurePrime
    (
        alpha,
        alphaMax_
    );

    nuFric_ = dimensionedScalar("0", nuFric_.dimensions(), 0.0);
    forAll(phaseIndexes_, phasei)
    {
        const phaseModel& phase = fluid_.phases()[phaseIndexes_[phasei]];
        tmp<volTensorField> tgradU(fvc::grad(phase.U()));
        const volTensorField& gradU(tgradU());
        volSymmTensorField D(symm(gradU));

        nuFric_ += frictionalStressModel_->nu
        (
            phase,
            alpha,
            alphaMax_,
            phase*Pfr_/phase.rho(),
            D
        )*phase;
    }
    nuFric_ /= max(alpha, residualAlpha_);
}


void Foam::kineticTheorySystem::correctAlphap()
{
    if (!alphapPtr_.valid())
    {
        return;
    }

    alphapPtr_() = 0.0;
    forAll(packingPhases_, phasei)
    {
        alphapPtr_() += packingPhases_[phasei];
    }
    alphapPtr_().correctBoundaryConditions();
}
// ************************************************************************* //
