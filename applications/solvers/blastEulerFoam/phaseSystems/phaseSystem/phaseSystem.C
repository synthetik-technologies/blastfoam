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
#include "interfacialPressureModel.H"
#include "interfacialVelocityModel.H"
#include "dragModel.H"
#include "dragODE.H"
#include "surfaceInterpolate.H"
#include "fvcDdt.H"

#include "SortableList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseSystem, 0);
}
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
                interfacialVelocityIter()->Ui()
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
                    if (!phase1.granular())
                    {
                        phase1.alphaRhoE() += deltaM & Ui;
                    }
                    phase2.alphaRhoU(nodej) -= deltaM;
                    if (!phase2.granular())
                    {
                        phase2.alphaRhoE() -= deltaM & Ui;
                    }
                }
            }
        }

        if
        (
            (phase1.granular() && phase1.nNodes() == 1 && !phase2.granular())
         || (!phase1.granular() && phase2.granular() && phase2.nNodes() == 1)
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
                if (!phase1.granular())
                {
                    phase1.alphaRhoE() += Fl & Ui;
                }
                phase2.alphaRhoU(nodej) -= Fl;
                if (!phase2.granular())
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
                      + (phase1.U() & fvc::grad(phase1.U()))
                      - fvc::ddt(phase2.U())
                      - (phase2.U() & fvc::grad(phase2.U()))
                    )
                );

                const volVectorField& Ui(*Uis[pair]);
                phase1.alphaRhoU(nodei) += Fvm;
                if (!phase1.granular())
                {
                    phase1.alphaRhoE() += Fvm & Ui;
                }
                phase2.alphaRhoU(nodej) -= Fvm;
                if (!phase2.granular())
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
                if (!phase1.granular())
                {
                    phase1.alphaRhoE() += Fwl & Ui;
                }
                phase2.alphaRhoU(nodej) -= Fwl;
                if (!phase2.granular())
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
                if (!phase1.granular())
                {
                    phase1.alphaRhoE() += Fwl & Ui;
                }
                phase2.alphaRhoU(nodej) -= Fwl;
                if (!phase2.granular())
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


void Foam::phaseSystem::calcMixtureVariables()
{
    tmp<volVectorField> alphaRhoU;
    tmp<volScalarField> alphaRhop;
    tmp<volScalarField> alphaRhoT;
    tmp<volScalarField> palphaRho;
    {
        const phaseModel& phase = phaseModels_[0];
        rho_ = phase.alphaRho();
        alphaRhoU = phase.alphaRho()*phase.U();
        alphaRhoT = phase.alphaRho()*phase.T();
        if (!phase.granular())
        {
            palphaRho =
                tmp<volScalarField>(new volScalarField(phase.alphaRho()));
            alphaRhop = phase.alphaRho()*phase.p();
        }
    }
    for (label phasei = 1; phasei < phaseModels_.size(); ++ phasei)
    {
        const phaseModel& phase = phaseModels_[phasei];
        const volScalarField& alphaRho = phase.alphaRho();
        rho_ += alphaRho;
        alphaRhoU.ref() += alphaRho*phase.U();
        alphaRhoT.ref() += alphaRho*phase.T();
        if (!phase.granular())
        {
            if (alphaRhop.valid())
            {
                palphaRho.ref() += alphaRho;
                alphaRhop.ref() += alphaRho*phase.p();
            }
            else
            {
                palphaRho =
                    tmp<volScalarField>(new volScalarField(alphaRho));
                alphaRhop = alphaRho*phase.p();
            }
        }
    }

    p_ = alphaRhop/palphaRho;
    U_ = alphaRhoU/rho_;
    T_ = alphaRhoT/rho_;
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
    integrationSystem("fluid", mesh),

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
        dimensionedVector("0", dimVelocity, Zero)
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

    g_(mesh.lookupObject<uniformDimensionedVectorField>("g")),

    phaseModels_(lookup("phases"), phaseModel::iNew(*this)),

    kineticTheoryPtr_(NULL),
    polydisperseKineticTheory_(false)
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
    generatePairsAndSubModels("aspectRatio", aspectRatioModels_);
    generatePairsAndSubModels("drag", dragModels_);
    generatePairsAndSubModels("virtualMass", virtualMassModels_);
    generatePairsAndSubModels("lift", liftModels_);
    generatePairsAndSubModels("turbulentDispersion", turbulentDispersionModels_);
    generatePairsAndSubModels("wallLubrication", wallLubricationModels_);
    generatePairsAndSubModels("heatTransfer", heatTransferModels_);
    generatePairsAndSubModels
    (
        "interfacialPressure",
        interfacialPressureModels_
    );
    generatePairsAndSubModels
    (
        "interfacialVelocity",
        interfacialVelocityModels_
    );

    if (lookupOrDefault<Switch>("solveDragODE", false))
    {
        dragODE_.set(new dragODE(*this, dragModels_));
    }

    if (phaseModels_.size() == 2)
    {
        phaseModels_.last().solveAlpha(false);
    }
    else
    {
        forAll(phaseModels_, phasei)
        {
            phaseModels_[phasei].solveAlpha(true);
        }
    }

    volScalarField sumAlpha("sumAlpha", phaseModels_[0]);
    label nFluids = phaseModels_[0].granular() ? 0 : 1;
    for (label phasei = 1; phasei < phaseModels_.size(); phasei++)
    {
        // Update boundaries
        phaseModels_[phasei].correctBoundaryConditions();

        sumAlpha += phaseModels_[phasei];
        if (!phaseModels_[phasei].granular())
        {
            nFluids++;
        }
    }

    if (nFluids > 1)
    {
        FatalErrorInFunction
            << "Only one fluid phase is currently supported by blastEulerFoam. "
            << "Multifluid implementations are under development."
            << endl
            << abort(FatalError);
    }

    if
    (
        max(phaseModels_.last()).value() == 0
     && min(phaseModels_.last()).value() == 0
    )
    {
        refCast<volScalarField>(phaseModels_.last()) = 1.0 - sumAlpha;
    }
    else if
    (
        max(sumAlpha()).value() - 1 > small
     && min(sumAlpha()).value() - 1 > small
    )
    {
        Info<<max(sumAlpha()).value()<<' '<<min(sumAlpha()).value()<<endl;
        FatalErrorInFunction
            << "Initial volume fractions do not sum to one." << nl
            << "min(sum(alphas)) = " << min(sumAlpha).value()
            << ", max(sum(alphas)) = " << max(sumAlpha).value()
            << endl
            << abort(FatalError);
    }

    encode();

    // Check if a granular phase is used and store at pointer if it is
    if (mesh_.foundObject<kineticTheorySystem>("kineticTheorySystem"))
    {
        kineticTheoryPtr_ =
        (
            &mesh_.lookupObjectRef<kineticTheorySystem>("kineticTheorySystem")
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
        volScalarField& alpha2(phaseModels_[1]);
        phaseModels_[0].decode();
        alpha2 = 1.0 - alpha1;
        alpha2.correctBoundaryConditions();
        phaseModels_[1].decode();
    }
    else
    {
        // find largest volume fraction and set to 1-sum
        forAll(rho_, celli)
        {
            SortableList<scalar> alphas(phaseModels_.size());
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
                alphas[phasei] = phaseModels_[phasei][celli];
            }
            alphas.reverseSort();

            const label fixedPhase = alphas.indices()[0];

            scalar sumAlpha = 0.0;
            for (label phasei = 1; phasei < alphas.size(); phasei++)
            {
                sumAlpha += alphas[phasei];
            }
            phaseModels_[fixedPhase][celli] = 1.0 - sumAlpha;
        }

        forAll(phaseModels_, phasei)
        {
            phaseModels_[phasei].correctBoundaryConditions();
            phaseModels_[phasei].decode();
        }
    }

    if (kineticTheoryPtr_ != NULL)
    {
        kineticTheoryPtr_->correct();
    }

    calcMixtureVariables();
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
        phaseModels_[phasei].postUpdate();
    }

    const dimensionedScalar& deltaT(mesh_.time().deltaT());
    relaxVelocity(deltaT);
    relaxTemperature(deltaT);
}


void Foam::phaseSystem::clearODEFields()
{
    Info<< "Phase statistics:"<<endl;
    forAll(phaseModels_, phasei)
    {
        Info<< phaseModels_[phasei].name() << ":" << endl;
        phaseModels_[phasei].clearODEFields();

        const volScalarField& alpha(phaseModels_[phasei]);
        const volScalarField& T(phaseModels_[phasei].T());
        Info<< "\t"
            << alpha.name() << " fraction, min, max = "
            << alpha.weightedAverage(mesh_.V()).value()
            << ' ' << min(alpha).value()
            << ' ' << max(alpha).value() << nl
            << "\t"
            << T.name() << " min, max = "
            << ' ' << min(T).value()
            << ' ' << max(T).value()
            << endl;
        if (!phaseModels_[phasei].granular())
        {
            const volScalarField& p(phaseModels_[phasei].p());
            Info<< "\t"
                << p.name() << " min, max = "
                << ' ' << min(p).value()
                << ' ' << max(p).value() << endl;
        }
        Info<< endl;

    }
    if (kineticTheoryPtr_)
    {
        if (kineticTheoryPtr_->polydisperse())
        {
            const volScalarField& alpha(kineticTheoryPtr_->alphap());
            Info<< alpha.name() << " fraction, min, max = "
                << alpha.weightedAverage(mesh_.V()).value()
                << ' ' << min(alpha).value()
                << ' ' << max(alpha).value()
                << endl;
        }
    }
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
        return aspectRatioModels_[key]->E(celli, nodei, nodej);
    }
    else
    {
        return 1.0;
    }
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
