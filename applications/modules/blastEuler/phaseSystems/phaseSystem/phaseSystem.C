/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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

#include "phaseSystem.H"
#include "kineticTheorySystem.H"
#include "aspectRatioModel.H"
#include "dragModel.H"
#include "virtualMassModel.H"
#include "wallLubricationModel.H"
#include "liftModel.H"
#include "turbulentDispersionModel.H"
#include "heatTransferModel.H"
#include "massTransferModel.H"
#include "interfacialPressureModel.H"
#include "interfacialVelocityModel.H"
#include "pressureRelaxationModel.H"
#include "dragModel.H"
#include "dragODE.H"
#include "pressureRelaxationSolver.H"
#include "surfaceInterpolate.H"
#include "fvcDdt.H"
#include "phaseFluxScheme.H"
#include "multicomponentBlastThermo.H"

#include "SortableList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseSystem, 0);
}

const Foam::dimensionedScalar
Foam::phaseSystem::zeroMDot(dimDensity/dimTime, 0.0);

template<>
const char* Foam::NamedEnum<Foam::phaseSystem::PIPressure, 2>::names[] =
    {
        "volume",
        "total"
    };

template<>
const char* Foam::NamedEnum<Foam::phaseSystem::PVRelaxation, 4>::names[] =
    {
        "none",
        "model",
        "ode",
        "instant"
    };

const Foam::NamedEnum<Foam::phaseSystem::PIPressure, 2>
Foam::phaseSystem::PIPressureNames_;

const Foam::NamedEnum<Foam::phaseSystem::PVRelaxation, 4>
Foam::phaseSystem::PVRelaxationNames_;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::phaseSystem::generatePairs
(
    const dictTable& modelDicts
)
{
    forAllConstIter(dictTable, modelDicts, iter)
    {
        const phasePairKey& key = iter.key();

        // pair already exists
        if (phasePairs_.found(key))
        {}

        // new ordered pair
        else if (key.ordered())
        {
            phasePairs_.insert
            (
                key,
                autoPtr<phasePair>
                (
                    new orderedPhasePair
                    (
                        phaseModels_[key.first()],
                        phaseModels_[key.second()]
                    )
                )
            );
        }

        // new unordered pair
        else
        {
            phasePairs_.insert
            (
                key,
                autoPtr<phasePair>
                (
                    new phasePair
                    (
                        phaseModels_[key.first()],
                        phaseModels_[key.second()]
                    )
                )
            );
        }
    }
}


#define StabAlphaRho(Phase, var)                        \
if (!stabAlphaRhos.set(Phase.index()))                  \
{                                                       \
    stabAlphaRhos.set                                   \
    (                                                   \
        Phase.index(),                                  \
        max(Phase.alphaRho(), Phase.residualAlphaRho()) \
    );                                                  \
}                                                       \
const volScalarField& var = stabAlphaRhos[Phase.index()];

void Foam::phaseSystem::relaxVelocity(const dimensionedScalar& deltaT)
{
    if (VRelaxation_ == NONE)
    {
        return;
    }

    if (dragODE_.valid())
    {
        Info<< "Solving drag ODE system" <<endl;
        dragODE_->solve(deltaT.value());
    }
    else if (VRelaxation_ == INSTANT)
    {
        volVectorField VI(fluidPhaseModels_[0].alphaRhoU());
        volScalarField fluidRho(fluidPhaseModels_[0].alphaRho());
        for (label i = 1; i < fluidPhaseModels_.size(); i++)
        {
            VI += fluidPhaseModels_[i].alphaRhoU();
            fluidRho += fluidPhaseModels_[i].alphaRho();
        }
        fluidRho.max(1e-6);
        VI /= fluidRho;

        forAll(fluidPhaseModels_, i)
        {
            phaseModel& phase = fluidPhaseModels_[i];
            phase.alphaRhoE() +=
                0.5*phase.alphaRho()*(magSqr(VI) - magSqr(phase.U()));
            phase.alphaRhoU() = phase.alphaRho()*VI;// - phase.U());
            phase.U() = VI;
        }
    }

    UiTable Uis;
    forAllConstIter
    (
        interfacialVelocityModelTable,
        interfacialVelocityModels_,
        interfacialVelocityIter
    )
    {
        const phasePair& pair(this->phasePairs_[interfacialVelocityIter.key()]);

        Uis.insert
        (
            pair,
            new volVectorField
            (
                IOobject::groupName("Ui", pair.name()),
                interfacialVelocityIter()->UI()
            )
        );
    }
    PtrList<volScalarField> stabAlphaRhos(phaseModels_.size());

    forAllConstIter
    (
        dragModelTable,
        dragModels_,
        dragModelIter
    )
    {
        const phasePair& pair(this->phasePairs_[dragModelIter.key()]);
        phaseModel& phase1 = phaseModels_[pair.phase1().name()];
        phaseModel& phase2 = phaseModels_[pair.phase2().name()];

        volScalarField Kd(dragModelIter()->K());

        StabAlphaRho(phase1, alphaRho1);

        if (!dragODE_.valid())
        {
            StabAlphaRho(phase2, alphaRho2);

            // Momentum and heat transfer
            volScalarField XiD(1.0/alphaRho1 + 1.0/alphaRho2);

            volVectorField deltaM
            (
                (phase1.U() - phase2.U())/XiD
                *(1.0/(Kd*XiD*deltaT + 1.0) - 1.0)
            );

            phase1.alphaRhoU() += deltaM;
            if (phase1.totalEnergy())
            {
                phase1.alphaRhoE() += deltaM & (*Uis[pair]);
            }

            phase2.alphaRhoU() -= deltaM;
            if (phase2.totalEnergy())
            {
                phase2.alphaRhoE() -= deltaM & (*Uis[pair]);
            }
        }


        if
        (
            (phase1.granular() && !phase2.granular())
         || (!phase1.granular() && phase2.granular())
        )
        {
            phaseModel* particles;
            phaseModel* gas;
            if (phase1.granular())
            {
                particles = &phase1;
                gas = &phase2;
            }
            else
            {
                particles = &phase2;
                gas = &phase1;
            }
            StabAlphaRho((*particles), alphaRhop);
            StabAlphaRho((*gas), alphaRhog);

            volScalarField XiD
            (
                1.0/alphaRhop + 1.0/alphaRhog
            );
            volScalarField ThetaOld
            (
                particles->alphaRhoPTE()/(1.5*alphaRhop)
            );
            volScalarField ThetaStar
            (
                ThetaOld*exp(-2.0*Kd*deltaT/alphaRhop)
            );
            ThetaStar.max(0.0);
            volScalarField ThetaStarStar
            (
                pow
                (

                    particles->productionSource(*gas)
                   /(Kd*XiD*deltaT + 1.0)
                   /alphaRhop
                   *deltaT
                  + pow(ThetaStar, 1.5),
                    2.0/3.0
                )
            );


            gas->alphaRhoE() -=
                1.5*particles->alphaRho()*(ThetaStarStar - ThetaOld);
            particles->alphaRhoPTE() =
                1.5*particles->alphaRho()*ThetaStarStar;
        }
    }

    //- Transformation of granular energy to thermal energy
    //  due to inelastic collisions
    forAll(phaseModels_, phasei)
    {
        phaseModel& phase1 = phaseModels_[phasei];
        if (phase1.granular())
        {
            forAll(phaseModels_, phasej)
            {
                phaseModel& phase2 = phaseModels_[phasej];
                if (phase2.granular())
                {
                    volScalarField gammaDot
                    (
                        phase1.dissipationSource
                        (
                            phase2,
                            mesh_.time().deltaT()
                        )
                    );

                    phase1.alphaRhoPTE() += gammaDot;
                    phase1.alphaRhoE() -= gammaDot;

                    if (phasei != phasej)
                    {
                        phase2.alphaRhoPTE() += gammaDot;
                        phase2.alphaRhoE() -= gammaDot;
                    }
                }
            }
        }
    }

    forAllConstIter
    (
        liftModelTable,
        liftModels_,
        liftModelIter
    )
    {
        const phasePair& pair(this->phasePairs_[liftModelIter.key()]);
        phaseModel& phase1 = phaseModels_[pair.phase1().name()];
        phaseModel& phase2 = phaseModels_[pair.phase2().name()];

        volVectorField Fl(liftModelIter()->F()*deltaT);

        phase1.alphaRhoU() += Fl;
        if (phase1.totalEnergy())
        {
            phase1.alphaRhoE() += Fl & (*Uis[pair]);
        }
        phase2.alphaRhoU() -= Fl;
        if (phase2.totalEnergy())
        {
            phase2.alphaRhoE() -= Fl & (*Uis[pair]);
        }
    }

    forAllConstIter
    (
        virtualMassModelTable,
        virtualMassModels_,
        virtualMassIter
    )
    {
        const phasePair& pair(this->phasePairs_[virtualMassIter.key()]);
        phaseModel& phase1 = phaseModels_[pair.phase1().name()];
        phaseModel& phase2 = phaseModels_[pair.phase2().name()];

        volVectorField Fvm
        (
          - virtualMassIter()->K()*deltaT
           *(
                fvc::ddt(phase1.U())
              + fvc::div(phase1.phi(), phase1.U())
              - fvc::ddt(phase2.U())
              - fvc::div(phase2.phi(), phase2.U())
            )
        );

        phase1.alphaRhoU() += Fvm;
        if (phase1.totalEnergy())
        {
            phase1.alphaRhoE() += Fvm & (*Uis[pair]);
        }
        phase2.alphaRhoU() -= Fvm;
        if (phase2.totalEnergy())
        {
            phase2.alphaRhoE() -= Fvm & (*Uis[pair]);
        }
    }

    forAllConstIter
    (
        wallLubricationModelTable,
        wallLubricationModels_,
        wallLubricationIter
    )
    {
        const phasePair& pair(this->phasePairs_[wallLubricationIter.key()]);
        phaseModel& phase1 = phaseModels_[pair.phase1().name()];
        phaseModel& phase2 = phaseModels_[pair.phase2().name()];

        volVectorField Fwl
        (
            wallLubricationIter()->F()*deltaT
        );

        phase1.alphaRhoU() += Fwl;
        if (phase1.totalEnergy())
        {
            phase1.alphaRhoE() += Fwl & (*Uis[pair]);
        }
        phase2.alphaRhoU() -= Fwl;
        if (phase2.totalEnergy())
        {
            phase2.alphaRhoE() -= Fwl & (*Uis[pair]);
        }
    }

    forAllConstIter
    (
        turbulentDispersionModelTable,
        turbulentDispersionModels_,
        turbulentDispersionIter
    )
    {
        const phasePair& pair(this->phasePairs_[turbulentDispersionIter.key()]);
        phaseModel& phase1 = phaseModels_[pair.phase1().name()];
        phaseModel& phase2 = phaseModels_[pair.phase2().name()];

        volVectorField Fwl
        (
            turbulentDispersionIter()->D()*deltaT
            *phase1.gradAlpha()

        );

        phase1.alphaRhoU() += Fwl;
        if (phase1.totalEnergy())
        {
            phase1.alphaRhoE() += Fwl & (*Uis[pair]);
        }
        phase2.alphaRhoU() -= Fwl;
        if (phase2.totalEnergy())
        {
            phase2.alphaRhoE() -= Fwl & (*Uis[pair]);
        }
    }
}


void Foam::phaseSystem::relaxTemperature(const dimensionedScalar& deltaT)
{
    PtrList<volScalarField> stabAlphaRhos(phaseModels_.size());

    // Update thermal energy due to heat transfer
    forAllConstIter
    (
        heatTransferModelTable,
        heatTransferModels_,
        heatTransferIter
    )
    {
        const phasePair& pair(this->phasePairs_[heatTransferIter.key()]);
        phaseModel& phase1 = phaseModels_[pair.phase1().name()];
        phaseModel& phase2 = phaseModels_[pair.phase2().name()];

        const blendedHeatTransferModel& ht = *heatTransferModels_[pair];
        volScalarField Kh(ht.K());

        StabAlphaRho(phase1, alphaRho1);
        StabAlphaRho(phase2, alphaRho2);

        volScalarField Xie
        (
            1.0/(alphaRho1*phase1.thermo().Cv())
          + 1.0/(alphaRho2*phase2.thermo().Cv())
        );

        volScalarField deltaE
        (
            (phase1.Ts() - phase2.Ts())/Xie
           *(exp(-Kh*Xie*deltaT) - 1.0)
        );

        phase1.alphaRhoE() += deltaE;
        phase2.alphaRhoE() -= deltaE;
    }
}


void Foam::phaseSystem::relaxPressure(const dimensionedScalar& deltaT)
{
    if (pressureSolver_.valid())
    {
        Info<< "Solving pressure relaxation" <<endl;
        if (pressureSolver_->solve(deltaT.value()))
        {
            decode();
        }
    }
}


void Foam::phaseSystem::calcMixtureVariables()
{
    rho_ = Zero;
    volVectorField alphaRhoU
    (
        volVectorField::New
        (
            "alphaRhoU",
            mesh_,
            dimensionedVector(dimDensity*dimVelocity, Zero)
        )
    );
    volScalarField alphaRhoT
    (
        volScalarField::New
        (
            "alphaRhoT",
            mesh_,
            dimensionedScalar(dimDensity*dimTemperature, 0.0)
        )
    );
    forAll(phaseModels_, phasei)
    {
        const phaseModel& phase = phaseModels_[phasei];
        const volScalarField& alphaRho = phase.alphaRho();
        rho_ += alphaRho;
        alphaRhoU += phase.alphaRhoU();
        alphaRhoT += alphaRho*phase.T();
    }
    U_ = alphaRhoU/rho_;
    T_ = alphaRhoT/rho_;

    if (fluidPhaseModels_.size() == 1)
    {
        p_ = fluidPhaseModels_[0].p();
        return;
    }

    volScalarField sumAlpha(volScalarField::New("sumAlpha", mesh_, 0));
    volScalarField sumAlphaRho
    (
        volScalarField::New
        (
            "sumAlphaRho",
            mesh_,
            dimensionedScalar("0", dimDensity, 0)
        )
    );
    volScalarField sumAlphaRhoKappa
    (
        volScalarField::New
        (
            "sumAlphaRhoKappa",
            mesh_,
            dimensionedScalar("0", dimDensity*kappa_.dimensions(), 0)
        )
    );
    p_ =  Zero;
    kappa_ = Zero;
    if (PIPtr_.valid())
    {
        PIPtr_() = Zero;
    }
    forAll(fluidPhaseModels_, phasei)
    {
        const phaseModel& phase(fluidPhaseModels_[phasei]);
        sumAlpha += phase;
        sumAlphaRho += phase.alphaRho();
        p_ += phase*phase.p();
        sumAlphaRhoKappa += phase.alphaRho()*phase.thermo().kappa();
        if (PIPtr_.valid())
        {
            PIPtr_() +=
                phase*phase.p()
              + phase.alphaRho()*magSqr(phase.U() - U_);
        }
    }
    p_ /= sumAlpha;
    kappa_ = sumAlphaRhoKappa/sumAlphaRho;
}


void Foam::phaseSystem::calcMixtureFluxes()
{
    rhoPhi_ = Zero;
    phi_ = Zero;
    forAll(phaseModels_, phasei)
    {
        const phaseModel& phase = phaseModels_[phasei];
        rhoPhi_ += phase.alphaRhoPhi();
        phi_ += phase.flux().alphaf()*phase.phi();
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseSystem::phaseSystem
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "phaseProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    timeIntegrationSystem("fluid", mesh),

    mesh_(mesh),

    rho_
    (
        IOobject
        (
            "rho",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimDensity, 0.0)
    ),

    U_
    (
        IOobject
        (
            "U",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimVelocity, Zero),
        "zeroGradient"
    ),

    phi_
    (
        IOobject
        (
            "phi",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::flux(U_)
    ),
    rhoPhi_
    (
        IOobject
        (
            "rhoPhi",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::flux(U_*rho_)
    ),

    p_
    (
        IOobject
        (
            "p",
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    PIModel_(VOLUME),

    T_
    (
        IOobject
        (
            "T",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimTemperature, 0.0)
    ),

    kappa_
    (
        IOobject
        (
            "kappa",
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedScalar
        (
            "kappa",
            dimPower/dimLength/dimTemperature,
            1.0
        )
    ),

    g_(mesh.lookupObject<uniformDimensionedVectorField>("g")),

    phaseModels_(lookup("phases"), phaseModel::iNew(*this)),

    master_(masterSystemList::New(mesh)),

    dragODE_(nullptr),

    VRelaxation_
    (
        this->found("VRelaxation")
      ? PVRelaxationNames_.read(this->lookup("VRelaxation"))
      : MODEL
    ),
    PRelaxation_(NONE)
{
    // Blending methods
    forAllConstIter(dictionary, subDict("blending"), iter)
    {
        blendingMethods_.insert
        (
            iter().keyword(),
            blendingMethod::New
            (
                iter().keyword(),
                iter().dict(),
                phaseModels_.toc()
            )
        );
    }

    // Sub-models
    generatePairsAndSubModels("aspectRatio", aspectRatioModels_, true);
    generateBlendedPairsAndSubModels("drag", dragModels_, true);
    generateBlendedPairsAndSubModels("virtualMass", virtualMassModels_, false);
    generateBlendedPairsAndSubModels("lift", liftModels_, false);
    generateBlendedPairsAndSubModels("turbulentDispersion", turbulentDispersionModels_, false);
    generateBlendedPairsAndSubModels("wallLubrication", wallLubricationModels_, false);
    generateBlendedPairsAndSubModels("heatTransfer", heatTransferModels_, true);
    generatePairsAndSubModels("massTransfer", massTransferModels_, false);
    generatePairsAndSubModels
    (
        "interfacialPressure",
        interfacialPressureModels_,
        true
    );
    generatePairsAndSubModels
    (
        "interfacialVelocity",
        interfacialVelocityModels_,
        true
    );

    generatePairsAndSubModels
    (
        "pressureRelaxation",
        pressureRelaxationModels_,
        false
    );

    label nFluids = 0;
    forAll(phaseModels_, phasei)
    {
        if (!phaseModels_[phasei].slavePressure())
        {
            fluidPhaseModels_.resize(nFluids + 1);
            fluidPhaseModels_.set
            (
                nFluids++,
                &phaseModels_[phasei]
            );
        }
        else
        {
            slavePhaseModels_.resize(slavePhaseModels_.size() + 1);
            slavePhaseModels_.set
            (
                slavePhaseModels_.size() - 1,
                &phaseModels_[phasei]
            );
        }
    }
    if (slavePhaseModels_.size())
    {
        if (!fluidPhaseModels_.size())
        {
            FatalErrorInFunction
                << "Only slave phase models are being used. "
                << "A fluid must also be used." << endl
                << abort(FatalError);
        }
        else if (&fluidPhaseModels_[0] == &phaseModels_[0])
        {
            FatalErrorInFunction
                << "A slave phase model should be the first phase." << nl
                << "Please switch the order of "
                << slavePhaseModels_[0].name() << " and "
                << fluidPhaseModels_[0].name() << endl
                << abort(FatalError);
        }
    }

    if (nFluids > 1)
    {
        if (this->found("PIModel"))
        {
            PIModel_ = PIPressureNames_.read(this->lookup("PIModel"));
        }
        if (PIModel_ == TOTAL)
        {
            PIPtr_.set
            (
                new volScalarField
                (
                    IOobject
                    (
                        "PI",
                        mesh.time().name(),
                        mesh
                    ),
                    mesh,
                    dimensionedScalar(dimPressure, 0.0)
                )
            );
        }

        pressureSolver_ = pressureRelaxationSolver::New
        (
            *this,
            interfacialPressureModels_,
            pressureRelaxationModels_
        );
        PRelaxation_ = PVRelaxationNames_[pressureSolver_->type()];

        if (VRelaxation_ == INSTANT)
        {
            DynamicList<phasePairKey> pairsToRemove(dragModels_.size());
            forAllConstIter(dragModelTable, dragModels_, iter)
            {
                if
                (
                    !iter()->phase1().slavePressure()
                 && !iter()->phase2().slavePressure()
                )
                {
                    pairsToRemove.append(iter.key());
                }
            }
            forAll(pairsToRemove, i)
            {
                dragModels_.erase(dragModels_.find(pairsToRemove[i]));
            }
        }
    }
    else
    {
        PRelaxation_ = NONE;
        if (VRelaxation_ == INSTANT)
        {
            VRelaxation_ = MODEL;
        }
    }

    if (VRelaxation_ == ODE)
    {
        dragODE_.set(new dragODE(*this, dragModels_));
    }

    bool limitInitialAlpha =
        this->lookupOrDefault("limitInitialVolumeFraction", true);
    if (phaseModels_.size() == 2)
    {
        scalar minAlpha1 = 0.0;
        scalar maxAlpha1 = 1.0;
        if (!phaseModels_[0].slavePressure())
        {
            minAlpha1 = phaseModels_[0].residualAlpha().value();
        }
        if (!phaseModels_[1].slavePressure())
        {
            maxAlpha1 = 1.0 - phaseModels_[1].residualAlpha().value();
        }

        if (phaseModels_[1].slavePressure())
        {
            phaseModels_[1].solveAlpha(true);
            phaseModels_[0].solveAlpha(false);

            if (limitInitialAlpha)
            {
                phaseModels_[1].maxMin(1.0 - maxAlpha1, 1.0 - minAlpha1);
            }
            dynamicCast<volScalarField>(phaseModels_[0]) ==
                1.0 - phaseModels_[1];
        }
        else
        {
            phaseModels_[0].solveAlpha(true);
            phaseModels_[1].solveAlpha(false);

            if (limitInitialAlpha)
            {
                phaseModels_[0].maxMin(minAlpha1, maxAlpha1);
            }
            dynamicCast<volScalarField>(phaseModels_[1]) ==
                1.0 - phaseModels_[0];
        }
    }
    else
    {
        volScalarField sumAlpha
        (
            volScalarField::New
            (
                "sumAlpha",
                mesh,
                0.0
            )
        );

        forAll(phaseModels_, phasei)
        {
            // Update boundaries
            phaseModels_[phasei].solveAlpha(true);
            phaseModels_[phasei].correctBoundaryConditions();
            sumAlpha += phaseModels_[phasei];
        }

        if
        (
            max(phaseModels_.last()).value() == 0
         && min(phaseModels_.last()).value() == 0
        )
        {
            dynamicCast<volScalarField>(phaseModels_.last()) =
                1.0 - sumAlpha;
        }
        else if
        (
            max(sumAlpha()).value() - 1 > small
         && min(sumAlpha()).value() - 1 > small
        )
        {
            FatalErrorInFunction
                << "Initial volume fractions do not sum to one." << nl
                << "min(sum(alphas)) = " << min(sumAlpha).value()
                << ", max(sum(alphas)) = " << max(sumAlpha).value()
                << endl
                << "Maximum deviation from unity: "
                << max(mag(sumAlpha - 1.0)).value() << endl
                << abort(FatalError);
        }
    }

    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].initializeModels();
    }
    encode();

    // Initialize master systems
    master_.initialize();

    hasMassTransfer_.setSize
    (
        phaseModels_.size(),
        boolList(phaseModels_.size(), false)
    );
    forAllConstIter
    (
        massTransferModelTable,
        massTransferModels_,
        massTransferIter
    )
    {
        const phasePair& pair(this->phasePairs_[massTransferIter.key()]);
        hasMassTransfer_[pair.phase1().index()][pair.phase1().index()] = true;
        hasMassTransfer_[pair.phase2().index()][pair.phase2().index()] = true;
        hasMassTransfer_[pair.phase1().index()][pair.phase2().index()] = true;
    	hasMassTransfer_[pair.phase1().index()][pair.phase2().index()] = true;

        if (!pair.ordered())
        {
            orderedPhasePair key1(pair.phase1(), pair.phase2());
            mDots_.insert
            (
                key1,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName("mDot", key1.name()),
                        this->mesh().time().name(),
                        this->mesh()
                    ),
                    this->mesh(),
                    dimensionedScalar(dimDensity/dimTime, 0)
                )
            );
            // orderedPhasePair key2(pair.phase2(), pair.phase1());
            // mDots_.insert
            // (
            //     key2,
            //     new volScalarField
            //     (
            //         IOobject
            //         (
            //             IOobject::groupName("mDot", key2.name()),
            //             this->mesh().time().name(),
            //             this->mesh()
            //         ),
            //         this->mesh(),
            //         dimensionedScalar(dimDensity/dimTime, 0)
            //     )
            // );
        }
        else
        {
            mDots_.insert
            (
                pair,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName("mDot", pair.name()),
                        this->mesh().time().name(),
                        this->mesh()
                    ),
                    this->mesh(),
                    dimensionedScalar(dimDensity/dimTime, 0)
                )
            );
        }
    }

    calcMixtureVariables();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseSystem::~phaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseSystem::decode()
{
    if (phaseModels_.size() == 2)
    {
        phaseModel& phase1 = phaseModels_[0];
        phaseModel& phase2 = phaseModels_[1];

        if (phase2.slavePressure())
        {
            forAll(phase2, celli)
            {
                phase1.correctVolumeFraction(1.0 - phase2[celli], celli);
            }
        }
        else
        {
            phase1.correctBoundaryConditions();
            forAll(phase1, celli)
            {
                phase2.correctVolumeFraction(1.0 - phase1[celli], celli);
            }
        }
    }
    else
    {
        forAll(rho_, celli)
        {
            scalar fixedAlpha = 0.0;
            forAll(slavePhaseModels_, i)
            {
                fixedAlpha += slavePhaseModels_[i][celli];
            }
            if (fixedAlpha >= 1)
            {
                forAll(slavePhaseModels_, i)
                {
                    slavePhaseModels_[i].scaleVolumeFraction
                    (
                        fixedAlpha,
                        celli
                    );
                }
                forAll(fluidPhaseModels_, i)
                {
                    fluidPhaseModels_[i].correctVolumeFraction
                    (
                        0.0,
                        celli
                    );
                }
            }
            else if (fluidPhaseModels_.size() == 1)
            {
                fluidPhaseModels_[0].correctVolumeFraction
                (
                    1.0 - fixedAlpha,
                    celli
                );
            }
            else
            {
                scalar fluidAlpha = 0.0;
                forAll(fluidPhaseModels_, i)
                {
                    fluidAlpha += max(fluidPhaseModels_[i][celli], 0.0);
                }
                scalar scale =
                    max((1.0 - fixedAlpha)/max(fluidAlpha, 1e-6), 0.0);
                if (scale > small)
                {
                    forAll(fluidPhaseModels_, i)
                    {
                        fluidPhaseModels_[i].scaleVolumeFraction
                        (
                            scale,
                            celli
                        );
                    }
                }
            }
        }
    }

    // Decode now that volume fraction has been calculated
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctBoundaryConditions();
        phaseModels_[phasei].decode();
    }

    // Update all masterSystems
    master_.update();

    // Update total quantities
    calcMixtureVariables();

    // if (PRelaxation_ == INSTANT)
    // {
    //     forAll(fluidPhaseModels_, i)
    //     {
    //         fluidPhaseModels_[i].p() = p();
    //     }
    // }
}


Foam::tmp<Foam::surfaceScalarField> Foam::phaseSystem::phi() const
{
    return phi_;
}


void Foam::phaseSystem::encode()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].encode();
    }
}


void Foam::phaseSystem::update()
{
    decode();
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].update();
    }

    calcMixtureFluxes();

    //- Update mass transfer rates
    forAllConstIter
    (
        massTransferModelTable,
        massTransferModels_,
        massTransferIter
    )
    {
        const phasePairKey& key = massTransferIter.key();
        if (!key.ordered())
        {
            phasePairKey key1(key.first(), key.second(), true);
            *mDots_[key1] = massTransferIter()->K();
        }
        else
        {
            *mDots_[key] = massTransferIter()->K();
        }
    }

    forAllIter
    (
        massTransferModelTable,
        massTransferModels_,
        massTransferIter
    )
    {
        const phaseModel& phase1 = massTransferIter()->pair().phase1();
        const phaseModel& phase2 = massTransferIter()->pair().phase2();

        blastThermo& thermo1 = phaseModels_[phase1.name()].thermo();
        blastThermo& thermo2 = phaseModels_[phase2.name()].thermo();

        const volScalarField& mDot = *mDots_[massTransferIter.key()];
        const List<word> species1(massTransferIter()->phase1Species());
        const List<word> species2(massTransferIter()->phase2Species());

        forAll(species1, i)
        {
            const word& specieName = species1[i];
            if (thermo1.containsSpecie(specieName))
            {
                dynamicCast<multicomponentBlastThermo>
                (
                    thermo1
                ).addDelta
                (
                    specieName,
                    massTransferIter()->Y(phase1, specieName)*mDot
                );
            }
        }

        forAll(species2, i)
        {
            const word& specieName = species2[i];
            if (thermo2.containsSpecie(specieName))
            {
                dynamicCast<multicomponentBlastThermo>
                (
                    thermo2
                ).addDelta
                (
                    specieName,
                    -massTransferIter()->Y(phase2, specieName)*mDot
                );
            }
        }
    }
}


void Foam::phaseSystem::solve()
{
    forAll(phaseModels_, phasei)
    {
        DebugInfo
            << "Solving " << phaseModels_[phasei].name() << ":" << endl;
        phaseModels_[phasei].solve();
    }

    master_.solve();

    if (VRelaxation_ == INSTANT)
    {
        volVectorField VI(fluidPhaseModels_[0].alphaRhoU());
        volScalarField fluidRho(fluidPhaseModels_[0].alphaRho());
        for (label i = 1; i < fluidPhaseModels_.size(); i++)
        {
            VI += fluidPhaseModels_[i].alphaRhoU();
            fluidRho += fluidPhaseModels_[i].alphaRho();
        }
        fluidRho.max(1e-6);
        VI /= fluidRho;

        forAll(fluidPhaseModels_, i)
        {
            phaseModel& phase = fluidPhaseModels_[i];
            volVectorField U(phase.alphaRhoU()/max(phase.alphaRho(), phase.residualAlphaRho()));
            phase.alphaRhoE() += 0.5*phase.alphaRho()*magSqr(VI - U);
            phase.alphaRhoU() = phase.alphaRho()*VI;
            phase.U() = VI;
        }
    }

    if (PRelaxation_ == INSTANT)
    {
        pressureSolver_->solve(time().deltaTValue());
    }

    if (PRelaxation_ == INSTANT || VRelaxation_ == INSTANT)
    {
        forAll(fluidPhaseModels_, i)
        {
            fluidPhaseModels_[i].correctDeltas();
        }
    }
}


void Foam::phaseSystem::postExplicit()
{
    decode();

    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].postExplicit();
    }

    const dimensionedScalar& deltaT(mesh_.time().deltaT());
    relaxVelocity(deltaT);
    relaxTemperature(deltaT);

    master_.postExplicit();
}


void Foam::phaseSystem::postImplicit()
{
    decode();

    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].postImplicit();
    }
    master_.postImplicit();

    relaxPressure(mesh_.time().deltaT());

    decode();
}


void Foam::phaseSystem::storeExplicit()
{
    // Clear flux schemes
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].storeExplicit();
    }
    master_.storeExplicit();
}


void Foam::phaseSystem::clear()
{
    // Clear flux schemes
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].clear();
    }
    master_.clear();
}


void Foam::phaseSystem::printInfo() const
{
    Info<< "Total statistics:" << endl << incrIndent
        << indent
        << T_.name() << " max, min = "
        << ' ' << max(T_).value()
        << ' ' << min(T_).value() << nl
        << indent
        << p_.name() << " max, min = "
        << ' ' << max(p_).value()
        << ' ' << min(p_).value()
        << endl << decrIndent;

    Info<< "Phase statistics:"<< endl << incrIndent;
    forAll(phaseModels_, phasei)
    {
        Info<< indent << phaseModels_[phasei].name() << ":"
            << endl << incrIndent;

        const volScalarField& alpha(phaseModels_[phasei]);
        const volScalarField& T(phaseModels_[phasei].T());
        Info<< indent
            << alpha.name() << " average, max, min = "
            << alpha.weightedAverage(mesh_.V()).value()
            << ' ' << max(alpha).value()
            << ' ' << min(alpha).value() << nl
            << indent
            << T.name() << " max, min = "
            << ' ' << max(T).value()
            << ' ' << min(T).value()
            << endl;
        if (!phaseModels_[phasei].slavePressure())
        {
            const volScalarField& p(phaseModels_[phasei].p());
            Info<< indent
                << p.name() << " max, min = "
                << ' ' << max(p).value()
                << ' ' << min(p).value() << endl;
        }

        tmp<volScalarField> tTs(phaseModels_[phasei].Ts());
        const volScalarField& Ts = tTs();
        if ((&T) != (&Ts))
        {
            Info<< indent
                << Ts.name() << " max, min = "
                << ' ' << max(Ts).value()
                << ' ' << min(Ts).value() << endl;
        }
        Info<< endl << decrIndent;
    }
    forAll(master_, i)
    {
        const masterSystem& sys = master_[i];
        if (sys.polydisperse())
        {
            const volScalarField& alpha = sys.alpha();
            const volScalarField& alphaMax = sys.alphaMax();
            tmp<volScalarField> talphaByAlphaMax(alpha/alphaMax);
            const volScalarField& alphaByAlphaMax = talphaByAlphaMax();

            Info<< indent << sys.name() << ":" << endl << incrIndent
                << indent << alpha.name() << " average, max, min = "
                << alpha.weightedAverage(mesh_.V()).value()
                << ' ' << max(alpha).value()
                << ' ' << min(alpha).value() << nl
                << indent << alphaMax.name() << " average, max, min = "
                << alphaMax.weightedAverage(mesh_.V()).value()
                << ' ' << max(alphaMax).value()
                << ' ' << min(alphaMax).value() << nl
                << indent << alpha.name() << " / " << alphaMax.name()
                << " average, max, min = "
                << alphaByAlphaMax.weightedAverage(mesh_.V()).value()
                << ' ' << max(alphaByAlphaMax).value()
                << ' ' << min(alphaByAlphaMax).value() << nl
                << decrIndent << endl;
        }
    }
    Info<< decrIndent;

}


Foam::scalar Foam::phaseSystem::cellPI(const label celli) const
{
    scalar PI = 0.0;
    if (PIModel_ == VOLUME)
    {
        forAll(fluidPhaseModels_, phasei)
        {
            const phaseModel& phase(fluidPhaseModels_[phasei]);
            PI += phase[celli]*phase.p()[celli];
        }
    }
    else if (PIModel_ == TOTAL)
    {
        vector VI = Zero;
        scalar rho = 0.0;
        forAll(fluidPhaseModels_, phasei)
        {
            const phaseModel& phase(fluidPhaseModels_[phasei]);
            VI += phase.alphaRhoU()[celli];
            rho += phase.alphaRho()[celli];
        }
        VI /= max(rho, 1e-6);

        forAll(fluidPhaseModels_, phasei)
        {
            const phaseModel& phase(fluidPhaseModels_[phasei]);
            PI +=
                phase[celli]*phase.p()[celli]
              + phase.alphaRho()[celli]
               *magSqr
                (
                    phase.alphaRhoU()[celli]
                   /max
                    (
                        phase.alphaRho()[celli],
                        phase.residualAlphaRho().value()
                    )
                  - VI
                );
        }

    }
    return PI;
}

Foam::tmp<Foam::volScalarField>
Foam::phaseSystem::E(const phasePairKey& key) const
{
    if (aspectRatioModels_.found(key))
    {
        return aspectRatioModels_[key]->E();
    }
    else
    {
        return volScalarField::New
        (
            aspectRatioModel::typeName + ":E",
            this->mesh_,
            dimensionedScalar(dimless, 1)
        );
    }
}


Foam::scalar Foam::phaseSystem::cellE
(
    const phasePairKey& key,
    const label celli
) const
{
    if (aspectRatioModels_.found(key))
    {
        return aspectRatioModels_[key]->cellE(celli);
    }
    else
    {
        return 1.0;
    }
}


bool Foam::phaseSystem::hasMassTransfer(const phaseModel& phase) const
{
    return hasMassTransfer_[phase.index()][phase.index()];
}


bool Foam::phaseSystem::hasMassTransfer
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    return
        phase1.index() != phase2.index()
     && (
            hasMassTransfer_[phase1.index()][phase2.index()]
         || hasMassTransfer_[phase2.index()][phase1.index()]
        );
}


Foam::tmp<Foam::volScalarField>
Foam::phaseSystem::mDot(const phaseModel& phase1, const phaseModel& phase2) const
{
    tmp<volScalarField> tmpmDoti
    (
        volScalarField::New
        (
            "mD" + phase1.name() + "." + phase2.name(),
            mesh_,
            dimensionedScalar(dimDensity/dimTime, 0.0)
        )
    );
    volScalarField::Internal& mDoti = tmpmDoti.ref();
    phasePairKey key1(phase1.name(), phase2.name(), true);
    phasePairKey key2(phase2.name(), phase1.name(), true);

    if (mDots_.found(key1))
    {
        mDoti += *mDots_[key1];
    }
    if (mDots_.found(key2))
    {
        mDoti -= *mDots_[key2];
    }
    return tmpmDoti;
}

Foam::tmp<Foam::volScalarField>
Foam::phaseSystem::mDotByRho
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    return mDotByRho(mDot(phase1, phase2)(), phase1, phase2);
}


Foam::tmp<Foam::volScalarField>
Foam::phaseSystem::mDotByRho
(
    const volScalarField& mD,
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    tmp<volScalarField> tmDotByRhoi
    (
        volScalarField::New
        (
            IOobject::groupName("mDotByRho", phase1.name()),
            mesh_,
            dimensionedScalar(inv(dimTime), 0.0)
        )
    );
    tmDotByRhoi.ref().internalFieldRef() =
        max(mD(), zeroMDot)/phase2.rho()()
      + min(mD(), zeroMDot)/phase1.rho()();
    return tmDotByRhoi;
}


Foam::tmp<Foam::volVectorField> Foam::phaseSystem::mDotU
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    return mDotU(mDot(phase1, phase2), phase1, phase2);
}


Foam::tmp<Foam::volVectorField> Foam::phaseSystem::mDotU
(
    const volScalarField& mD,
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    tmp<volVectorField> tmDotUi
    (
        volVectorField::New
        (
            IOobject::groupName("mDotU", phase1.name()),
            mesh_,
            dimensionedVector(dimDensity*dimVelocity/dimTime, Zero)
        )
    );
    tmDotUi.ref().internalFieldRef() =
        max(mD(), zeroMDot)*phase2.U()()
      + min(mD(), zeroMDot)*phase1.U()();
    return tmDotUi;
}

Foam::tmp<Foam::volScalarField> Foam::phaseSystem::mDotE
(
    const volScalarField& mD,
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    volScalarField::Internal mD21(max(mD(), zeroMDot));
    volScalarField::Internal mD12(min(mD(), zeroMDot));

    tmp<volScalarField> tmDotEi
    (
        volScalarField::New
        (
            IOobject::groupName("mDotE", phase1.name()),
            mesh_,
            dimensionedScalar(dimDensity*sqr(dimVelocity)/dimTime, 0.0)
        )
    );
    volScalarField::Internal& mDotEi = tmDotEi.ref();

    // Determine heat of formation of the reactants/products
    tmp<volScalarField> hc1, hc2;
    const basicThermo& bthermo1 = phase1.thermo();
    const basicThermo& bthermo2 = phase2.thermo();
    if
    (
        isA<multicomponentThermo>(bthermo1)
     || isA<multicomponentThermo>(bthermo2)
    )
    {
        phasePairKey key1(phase1.name(), phase2.name(), true);
        phasePairKey key2(phase2.name(), phase1.name(), true);

        const massTransferModel& mt =
            massTransferModels_.found(key1)
          ? massTransferModels_[key1]()
          : massTransferModels_[key2]();
        const List<word> species1(mt.species(phase1));
        const List<word> species2(mt.species(phase2));

        // Add heat of formation of the produced species
        if (isA<multicomponentThermo>(bthermo1))
        {
            const multicomponentThermo& thermo1 =
                dynamicCast<const multicomponentThermo>(bthermo1);
            forAll(species1, i)
            {
                const word& specieName = species1[i];
                const label speciei = thermo1.species()[specieName];
                if (hc1.valid())
                {
                    hc1.ref() += mt.Y(phase1, specieName)*thermo1.hfi(speciei);
                }
                else
                {
                    hc1 = mt.Y(phase1, specieName)*thermo1.hfi(speciei);
                }
            }
        }

        // Add heat of formation of the consumed species
        if (isA<multicomponentThermo>(bthermo2))
        {
            const multicomponentThermo& thermo2 =
                dynamicCast<const multicomponentThermo>(bthermo2);
            forAll(species2, i)
            {
                const word& specieName = species2[i];
                const label speciei = thermo2.species()[specieName];
                if (hc2.valid())
                {
                    hc2.ref() += mt.Y(phase2, specieName)*thermo2.hfi(speciei);
                }
                else
                {
                    hc2 = mt.Y(phase2, specieName)*thermo2.hfi(speciei);
                }
            }
        }
    }

    // No heat of formation set so use the mixture
    if (!hc1.valid())
    {
        hc1 = phase1.thermo().hc();
    }
    if (!hc2.valid())
    {
        hc2 = phase2.thermo().hc();
    }

    mDotEi =
        mD21*(phase2.he() + (hc2()() - hc1()()))
      + mD12*phase1.he()();

    // Add kinetic energy contributions
    if (phase1.totalEnergy())
    {
        tmp<volScalarField::Internal> K1(0.5*magSqr(phase1.U()()));
        tmp<volScalarField::Internal> K2(0.5*magSqr(phase2.U()()));

        if (phase2.granular())
        {
            K2.ref() += 1.5*phase2.Theta()();
        }

        mDotEi += mD21*K2 + mD12*K1;
    }

    return tmDotEi;
}


Foam::tmp<Foam::volScalarField>
Foam::phaseSystem::mDotE
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    return mDotE(mDot(phase1, phase2)(), phase1, phase2);
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::mDotPTE
(
    const volScalarField& mD,
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    tmp<volScalarField> tmDotPTEi
    (
        volScalarField::New
        (
            IOobject::groupName("mDotPTE", phase1.name()),
            mesh_,
            dimensionedScalar(dimDensity*sqr(dimVelocity)/dimTime, 0.0)
        )
    );
    volScalarField::Internal& mDotPTEi = tmDotPTEi.ref();

    if (phase1.granular())
    {
        mDotPTEi += min(mD(), zeroMDot)*1.5*phase1.Theta()();
    }
    if (phase2.granular())
    {
        mDotPTEi += max(mD(), zeroMDot)*1.5*phase2.Theta()();
    }
    return tmDotPTEi;
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::mDotPTE
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    return mDotPTE(mDot(phase1, phase2), phase1, phase2);
}



bool Foam::phaseSystem::read()
{
    IOdictionary& dict(*this);
    if (dict.regIOobject::read())
    {
        bool readOK = true;

        forAll(phaseModels_, phasei)
        {
            readOK &= phaseModels_[phasei].read();
        }

        // models ...

        return readOK;
    }
    return false;
}

// ************************************************************************* //
