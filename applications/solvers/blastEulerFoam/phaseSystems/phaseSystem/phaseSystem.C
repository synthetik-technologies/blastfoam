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


void Foam::phaseSystem::relaxVelocity(const dimensionedScalar& deltaT)
{
    if (dragODE_.valid())
    {
        Info<< "Solving drag ODE system" <<endl;
        dragODE_->solve(deltaT.value());
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

        volScalarField Kd
        (
            IOobject
            (
                IOobject::groupName("Kd", pair.name()),
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimDensity/dimTime, 0.0)
        );

        for (label nodei = 0; nodei < phase1.nNodes(); nodei++)
        {
            volScalarField alphaRho1
            (
                max(phase1.alphaRho(nodei), phase1.residualAlphaRho())
            );

            for (label nodej = 0; nodej < phase2.nNodes(); nodej++)
            {
                volScalarField dragCoeff
                (
                    dragModels_[pair]->K(nodei, nodej)
                );
                Kd += dragCoeff;

                if (!dragODE_.valid())
                {
                    volScalarField alphaRho2
                    (
                        max(phase2.alphaRho(nodej), phase2.residualAlphaRho())
                    );
                    // Momentum and heat transfer
                    volScalarField XiD
                    (
                        1.0/alphaRho1 + 1.0/alphaRho2
                    );

                    volVectorField deltaM
                    (
                        (phase1.U(nodei) - phase2.U(nodej))/XiD
                       *(1.0/(dragCoeff*XiD*deltaT + 1.0) - 1.0)
                    );

                    const volVectorField& Ui(*Uis[pair]);
                    phase1.alphaRhoU(nodei) += deltaM;
                    if (phase1.totalEnergy())
                    {
                        phase1.alphaRhoE() += deltaM & Ui;
                    }
                    phase2.alphaRhoU(nodej) -= deltaM;
                    if (phase2.totalEnergy())
                    {
                        phase2.alphaRhoE() -= deltaM & Ui;
                    }
                }
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
            volScalarField alphaRhop
            (
                max(particles->alphaRho(), particles->residualAlphaRho())
            );
            volScalarField alphaRhog
            (
                max(gas->alphaRho(), gas->residualAlphaRho())
            );
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
        forAll(phaseModels_, phasej)
        {
            if
            (
                phaseModels_[phasei].granular()
             && phaseModels_[phasej].granular()
            )
            {
                volScalarField gammaDot
                (
                    phaseModels_[phasei].dissipationSource
                    (
                        phaseModels_[phasej]
                    )*deltaT
                );
                phaseModels_[phasei].alphaRhoPTE() -= gammaDot;
                phaseModels_[phasei].alphaRhoE() += gammaDot;
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

        for (label nodei = 0; nodei < phase1.nNodes(); nodei++)
        {
            volScalarField alphaRho1
            (
                max(phase1.alphaRho(nodei), phase1.residualAlphaRho())
            );

            for (label nodej = 0; nodej < phase2.nNodes(); nodej++)
            {
                volScalarField alphaRho2
                (
                    max(phase2.alphaRho(nodej), phase2.residualAlphaRho())
                );

                volVectorField Fl
                (
                    liftModelIter()->F<vector>(nodei, nodej)*deltaT
                );

                const volVectorField& Ui(*Uis[pair]);
                phase1.alphaRhoU(nodei) += Fl;
                if (phase1.totalEnergy())
                {
                    phase1.alphaRhoE() += Fl & Ui;
                }
                phase2.alphaRhoU(nodej) -= Fl;
                if (phase2.totalEnergy())
                {
                    phase2.alphaRhoE() -= Fl & Ui;
                }
            }
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

        for (label nodei = 0; nodei < phase1.nNodes(); nodei++)
        {
            volScalarField alphaRho1
            (
                max(phase1.alphaRho(nodei), phase1.residualAlphaRho())
            );

            for (label nodej = 0; nodej < phase2.nNodes(); nodej++)
            {
                volScalarField alphaRho2
                (
                    max(phase2.alphaRho(nodej), phase2.residualAlphaRho())
                );

                volVectorField Fvm
                (
                   -virtualMassIter()->K(nodei, nodej)*deltaT
                   *(
                        fvc::ddt(phase1.U())
                      + fvc::div(phase1.phi(), phase1.U())
                      - fvc::ddt(phase2.U())
                      - fvc::div(phase2.phi(), phase2.U())
                    )
                );

                const volVectorField& Ui(*Uis[pair]);
                phase1.alphaRhoU(nodei) += Fvm;
                if (phase1.totalEnergy())
                {
                    phase1.alphaRhoE() += Fvm & Ui;
                }
                phase2.alphaRhoU(nodej) -= Fvm;
                if (phase2.totalEnergy())
                {
                    phase2.alphaRhoE() -= Fvm & Ui;
                }
            }
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

        for (label nodei = 0; nodei < phase1.nNodes(); nodei++)
        {
            volScalarField alphaRho1
            (
                max(phase1.alphaRho(nodei), phase1.residualAlphaRho())
            );

            for (label nodej = 0; nodej < phase2.nNodes(); nodej++)
            {
                volScalarField alphaRho2
                (
                    max(phase2.alphaRho(nodej), phase2.residualAlphaRho())
                );

                volVectorField Fwl
                (
                    wallLubricationIter()->F<vector>(nodei, nodej)*deltaT
                );

                const volVectorField& Ui(*Uis[pair]);
                phase1.alphaRhoU(nodei) += Fwl;
                if (phase1.totalEnergy())
                {
                    phase1.alphaRhoE() += Fwl & Ui;
                }
                phase2.alphaRhoU(nodej) -= Fwl;
                if (phase2.totalEnergy())
                {
                    phase2.alphaRhoE() -= Fwl & Ui;
                }
            }
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

        for (label nodei = 0; nodei < phase1.nNodes(); nodei++)
        {
            volScalarField alphaRho1
            (
                max(phase1.alphaRho(nodei), phase1.residualAlphaRho())
            );

            for (label nodej = 0; nodej < phase2.nNodes(); nodej++)
            {
                volScalarField alphaRho2
                (
                    max(phase2.alphaRho(nodej), phase2.residualAlphaRho())
                );

                volVectorField Fwl
                (
                    turbulentDispersionIter()->D(nodei, nodej)*deltaT
                   *phase1.gradAlpha()

                );

                const volVectorField& Ui(*Uis[pair]);
                phase1.alphaRhoU(nodei) += Fwl;
                if (phase1.totalEnergy())
                {
                    phase1.alphaRhoE() += Fwl & Ui;
                }
                phase2.alphaRhoU(nodej) -= Fwl;
                if (phase2.totalEnergy())
                {
                    phase2.alphaRhoE() -= Fwl & Ui;
                }
            }
        }
    }
}


void Foam::phaseSystem::relaxTemperature(const dimensionedScalar& deltaT)
{
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

        volScalarField KhMean
        (
            IOobject
            (
                "KhMean",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimensionSet(1, -1, -3, -1, 0), 0.0)
        );

        for (label nodei = 0; nodei < phase1.nNodes(); nodei++)
        {
            for (label nodej = 0; nodej < phase2.nNodes(); nodej++)
            {
                KhMean += heatTransferModels_[pair]->K(nodei, nodej);
            }
        }
        volScalarField alphaRho1
        (
            max(phase1.alphaRho(), phase1.residualAlphaRho())
        );
        volScalarField alphaRho2
        (
            max(phase2.alphaRho(), phase2.residualAlphaRho())
        );
        volScalarField Xie
        (
            1.0/(alphaRho1*phase1.Cv()) + 1.0/(alphaRho2*phase2.Cv())
        );
        volScalarField deltaE
        (
            (phase1.T() - phase2.T())/Xie
           *(exp(-KhMean*Xie*deltaT) - 1.0)
        );

        phase1.alphaRhoE() += deltaE;
        phase2.alphaRhoE() -= deltaE;
    }
}


void Foam::phaseSystem::relaxPressure(const dimensionedScalar& deltaT)
{
    if (pressureSolver_->solve())
    {
        Info<< "Solving pressure relaxation" <<endl;
        pressureSolver_->solve(deltaT.value());

        encode();
    }
}


void Foam::phaseSystem::calcMixtureVariables()
{
    rho_ = Zero;
    phi_ = Zero;
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
        alphaRhoU += alphaRho*phase.U();
        alphaRhoT += alphaRho*phase.T();
        phi_ += phase.alphaPhi();
    }
    U_ = alphaRhoU/rho_;
    T_ = alphaRhoT/rho_;

    if (fluidPhaseModels_.size() < 2)
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
        sumAlphaRhoKappa += phase.alphaRho()*phase.kappa();
        if (PIPtr_.valid())
        {
            PIPtr_() +=
                phase*(phase.p() + phase.rho()*magSqr(phase.U() - U_));
        }
    }
    p_ /= sumAlpha;
    kappa_ = sumAlphaRhoKappa/sumAlphaRho;
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
            IOobject::MUST_READ_IF_MODIFIED,
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
            mesh.time().timeName(),
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
            mesh.time().timeName(),
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
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::flux(U_)
    ),

    p_
    (
        IOobject
        (
            "p",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    T_
    (
        IOobject
        (
            "T",
            mesh.time().timeName(),
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
            mesh.time().timeName(),
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

    kineticTheoryPtr_(nullptr),
    polydisperseKineticTheory_(false),

    dragODE_(nullptr)
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
    generatePairsAndSubModels("drag", dragModels_, true);
    generatePairsAndSubModels("virtualMass", virtualMassModels_, false);
    generatePairsAndSubModels("lift", liftModels_, false);
    generatePairsAndSubModels("turbulentDispersion", turbulentDispersionModels_, false);
    generatePairsAndSubModels("wallLubrication", wallLubricationModels_, false);
    generatePairsAndSubModels("heatTransfer", heatTransferModels_, true);
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
        PIPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "PI",
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar(dimPressure, 0.0)
            )
        );
    }

    if (lookupOrDefault<Switch>("solveDragODE", false))
    {
        dragODE_.set(new dragODE(*this, dragModels_));
    }

    pressureSolver_ = pressureRelaxationSolver::New
    (
        *this,
        interfacialPressureModels_,
        pressureRelaxationModels_

    );

    if (phaseModels_.size() == 2)
    {
        if (phaseModels_[1].slavePressure())
        {
            phaseModels_[1].solveAlpha(true);
            phaseModels_[0].solveAlpha(false);
        }
        else
        {
            phaseModels_[0].solveAlpha(true);
            phaseModels_[1].solveAlpha(false);
        }
    }
    else
    {
        forAll(phaseModels_, phasei)
        {
            phaseModels_[phasei].solveAlpha(true);
        }
    }

    if (phaseModels_.size() == 2)
    {
        dynamicCast<volScalarField>(phaseModels_[1]) =
            1.0 - phaseModels_[0];
        phaseModels_[1].correctVolumeFraction();
    }
    else
    {
        volScalarField sumAlpha("sumAlpha", phaseModels_[0]);
        for (label phasei = 1; phasei < phaseModels_.size(); phasei++)
        {
            // Update boundaries
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

    // Check if a granular phase is used and store at pointer if it is
    if
    (
        mesh_.foundObject<kineticTheorySystem>
        (
            kineticTheorySystem::typeName
        )
    )
    {
        kineticTheoryPtr_.set
        (
            &mesh_.lookupObjectRef<kineticTheorySystem>
            (
                kineticTheorySystem::typeName
            )
        );

        // Initialize fields after all granular phases are initialized
        kineticTheoryPtr_->correct();

        //- If only one granular phase is used, the multiphase limiting is not
        //  needed so it is skipped
        if (kineticTheoryPtr_->polydisperse())
        {
            polydisperseKineticTheory_ = true;
        }
    }

    forAllConstIter
    (
        massTransferModelTable,
        massTransferModels_,
        massTransferIter
    )
    {
        const phasePair& pair(this->phasePairs_[massTransferIter.key()]);

        mDots_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("mDot", pair.name()),
                    this->mesh().time().timeName(),
                    this->mesh()
                ),
                this->mesh(),
                dimensionedScalar(dimDensity/dimTime, 0)
            )
        );
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
        volScalarField& alpha1(phaseModels_[0]);
        phaseModels_[0].correctVolumeFraction();

        volScalarField& alpha2(phaseModels_[1]);
        alpha2 = 1.0 - alpha1;
        phaseModels_[1].correctVolumeFraction();

        phaseModels_[0].decode();
        phaseModels_[1].decode();
    }
    else
    {
        forAll(phaseModels_, phasei)
        {
            phaseModels_[phasei].correctVolumeFraction();
        }

        // find largest volume fraction and set to 1-sum
        forAll(rho_, celli)
        {
            SortableList<scalar> alphas(phaseModels_.size(), 0.0);
            scalar sumAlpha = 0.0;
            forAll(phaseModels_, phasei)
            {
                phaseModels_[phasei][celli] =
                    Foam::max
                    (
                        Foam::min
                        (
                            phaseModels_[phasei][celli],
                            1.0
                        ),
                        0.0
                    );
                if (!phaseModels_[phasei].slavePressure())
                {
                    alphas[phasei] = phaseModels_[phasei][celli];
                }
                sumAlpha += phaseModels_[phasei][celli];
            }
            alphas.reverseSort();

            const label fixedPhase = alphas.indices()[0];
            phaseModels_[fixedPhase][celli] =
                1.0
              - (sumAlpha - phaseModels_[fixedPhase][celli]);
        }

        forAll(phaseModels_, phasei)
        {
            phaseModels_[phasei].correctVolumeFraction();
            phaseModels_[phasei].decode();
        }
    }

    if (kineticTheoryPtr_.valid())
    {
        kineticTheoryPtr_->correct();
    }

    calcMixtureVariables();

    forAll(phaseModels_, phasei)
    {
        phaseModel& phase(phaseModels_[phasei]);
        if (!phase.slavePressure())
        {
            phase.p() = p_;
        }
    }
    encode();
}


Foam::tmp<Foam::surfaceScalarField> Foam::phaseSystem::phi() const
{
    tmp<surfaceScalarField> phiTmp
    (
        new surfaceScalarField
        (
            "phi",
            phaseModels_[0].alphaRhoPhi()
        )
    );
    for (label phasei = 1; phasei < phaseModels_.size(); ++ phasei)
    {
        phiTmp.ref() += phaseModels_[phasei].alphaRhoPhi();
    }
    phiTmp.ref() /= fvc::interpolate(rho_);
    return phiTmp;
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

    //- Update mass transfer rates
    forAllConstIter
    (
        massTransferModelTable,
        massTransferModels_,
        massTransferIter
    )
    {
        *mDots_[massTransferIter.key()] = massTransferIter()->K();
    }

    forAllIter
    (
        massTransferModelTable,
        massTransferModels_,
        massTransferIter
    )
    {
        blastThermo& dispersedThermo
        (
            phaseModels_[massTransferIter()->pair().dispersed().name()].thermo()
        );
        blastThermo& continuousThermo
        (
            phaseModels_[massTransferIter()->pair().continuous().name()].thermo()
        );

        tmp<volScalarField> mDot(*mDots_[massTransferIter.key()]);
        List<word> species(massTransferIter()->dispersedSpecies());
        species.append(massTransferIter()->continuousSpecies());

        forAll(species, i)
        {
            const word& specieName(species[i]);
            if (isA<multicomponentBlastThermo>(dispersedThermo))
            {
                multicomponentBlastThermo& thermo =
                    dynamicCast<multicomponentBlastThermo>(dispersedThermo);
                if (thermo.contains(specieName))
                {
                    tmp<volScalarField> YmDot
                    (
                        -massTransferIter()->dispersedYi(specieName)*mDot
                    );
                    thermo.addDelta
                    (
                        specieName,
                        YmDot
                    );
                }
            }
            if (isA<multicomponentBlastThermo>(continuousThermo))
            {
                multicomponentBlastThermo& thermo =
                    dynamicCast<multicomponentBlastThermo>(continuousThermo);
                if (thermo.contains(specieName))
                {
                    tmp<volScalarField> YmDot
                    (
                        massTransferIter()->continuousYi(specieName)*mDot
                    );
                    thermo.addDelta
                    (
                        specieName,
                        YmDot
                    );
                }
            }
        }
    }
}


void Foam::phaseSystem::solve()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].solve();
    }
}


void Foam::phaseSystem::postUpdate()
{
    decode();
    forAll(phaseModels_, phasei)
    {
        Info<< "Solving: " << phaseModels_[phasei].name() << ":" << endl;
        phaseModels_[phasei].postUpdate();
        Info<< endl;
    }

    decode();

    const dimensionedScalar& deltaT(mesh_.time().deltaT());
    relaxVelocity(deltaT);
    relaxTemperature(deltaT);

    decode();

    relaxPressure(deltaT);
}


void Foam::phaseSystem::clear()
{
    // Clear flux schemes
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].flux().clear();
    }
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
        Info<< endl << decrIndent;
    }
    if (kineticTheoryPtr_.valid())
    {
        if (kineticTheoryPtr_->polydisperse())
        {
            const volScalarField& alpha(kineticTheoryPtr_->alphap());
            Info<< nl
                << indent << alpha.name() << " fraction, max, min = "
                << alpha.weightedAverage(mesh_.V()).value()
                << ' ' << max(alpha).value()
                << ' ' << min(alpha).value()
                << endl;
        }
    }
    Info<< decrIndent;

}


Foam::tmp<Foam::volScalarField>
Foam::phaseSystem::E
(
    const phasePairKey& key,
    const label nodei,
    const label nodej
) const
{
    if (aspectRatioModels_.found(key))
    {
        return aspectRatioModels_[key]->E(nodei, nodej);
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


Foam::scalar
Foam::phaseSystem::cellE
(
    const label celli,
    const phasePairKey& key,
    const label nodei,
    const label nodej
) const
{
    if (aspectRatioModels_.found(key))
    {
        return aspectRatioModels_[key]->cellE(celli, nodei, nodej);
    }
    else
    {
        return 1.0;
    }
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
    volScalarField& mDoti = tmpmDoti.ref();
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
    return mDotByRho(mDot(phase2, phase2)(), phase1, phase2);
}


Foam::tmp<Foam::volScalarField>
Foam::phaseSystem::mDotByRho
(
    const volScalarField& mD,
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    return
        volScalarField::New
        (
            IOobject::groupName("mDotByRho", phase1.name()),
            max(mD, zeroMDot)/phase2.rho() + min(mD, zeroMDot)/phase1.rho()
        );
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
    return
        volVectorField::New
        (
            IOobject::groupName("mDotU", phase1.name()),
            max(mD, zeroMDot)*phase2.U() + min(mD, zeroMDot)*phase1.U()
        );
}

Foam::tmp<Foam::volScalarField> Foam::phaseSystem::mDotE
(
    const volScalarField& mD,
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    tmp<volScalarField> tmDotEi
    (
        volScalarField::New
        (
            IOobject::groupName("mDotE", phase1.name()),
            mesh_,
            dimensionedScalar(dimDensity*sqr(dimVelocity)/dimTime, 0.0)
        )
    );
    volScalarField& mDotEi = tmDotEi.ref();

    volScalarField hc(phase1.thermo().hc() - phase2.thermo().hc());
    volScalarField mD21(max(mD, zeroMDot));
    volScalarField mD12(min(mD, zeroMDot));
    mDotEi += mD21*phase2.he() + mD12*phase1.he() + mD21*hc;

    if (phase1.totalEnergy())
    {
        volScalarField K1(0.5*magSqr(phase1.U()));
        volScalarField K2(0.5*magSqr(phase2.U()));

        if (phase2.granular())
        {
            K2 += 1.5*phase2.Theta();
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
    volScalarField& mDotPTEi = tmDotPTEi.ref();

    if (phase1.granular())
    {
        mDotPTEi += min(mD, zeroMDot)*1.5*phase1.Theta();
    }
    if (phase2.granular())
    {
        mDotPTEi += max(mD, zeroMDot)*1.5*phase2.Theta();
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
