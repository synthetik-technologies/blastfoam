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
#include "fluidBlastThermo.H"
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

void Foam::phaseSystem::relaxFluidVelocitiesInstant()
{
    const label nFluids = fluidPhaseModels_.size();
    UPtrList<fluidBlastThermo> thermos(nFluids);
    forAll(fluidPhaseModels_, i)
    {
        thermos.set
        (
            i,
            &mesh_.lookupObjectRef<fluidBlastThermo>
            (
                IOobject::groupName
                (
                    physicalProperties::typeName,
                    fluidPhaseModels_[i].name()
                )
            )
        );
    }

    scalarList m(nFluids), E(nFluids), Enew(nFluids);
    forAll(rho_, celli)
    {
        scalar M = 0;
        vector MU = Zero;
        forAll(fluidPhaseModels_, i)
        {
            const phaseModel& phase = fluidPhaseModels_[i];
            m[i] = phase.alphaRho()[celli];
            E[i] = phase.alphaRhoE()[celli];
            M += m[i];
            MU += phase.alphaRhoU()[celli];
        }
        const vector VI = MU/max(M, 1e-6);
        const scalar KI = 0.5*magSqr(VI);

        // Each phase keeps its internal energy and is heated by its own
        // dissipation, 0.5*alphaRho*|U - VI|^2, which conserves the total
        // energy. The kinetic energy of a phase with little mass can exceed
        // its total energy (its momentum comes from face fluxes weighted by
        // the face volume fraction); where the internal energy would fall
        // below the thermo range, that kinetic energy is not counted and the
        // small remainder goes to the other phases by mass. Heating beyond
        // the upper thermo limit is not applied.
        scalar R = 0;
        scalar Mvalid = 0;
        forAll(fluidPhaseModels_, i)
        {
            const phaseModel& phase = fluidPhaseModels_[i];
            if (m[i] <= phase.residualAlphaRho().value())
            {
                Enew[i] = E[i];
                continue;
            }
            const vector U = phase.alphaRhoU()[celli]/m[i];
            scalar e = (E[i] - 0.5*m[i]*magSqr(U))/m[i];
            scalar D = 0.5*magSqr(U - VI);
            const scalar eMin = thermos[i].cellhe(thermos[i].TLow(), celli);
            const scalar eMax = thermos[i].cellhe(thermos[i].THigh(), celli);
            if (e < eMin)
            {
                e = eMin;
                D = 0;
            }
            const scalar eNew = Foam::min(e + D, Foam::max(eMax, e));
            Enew[i] = m[i]*(eNew + KI);
            if (e + D <= eMax)
            {
                R += E[i] - Enew[i];
            }
            Mvalid += m[i];
        }

        forAll(fluidPhaseModels_, i)
        {
            phaseModel& phase = fluidPhaseModels_[i];
            scalar Ei = Enew[i];
            if (Mvalid > small && m[i] > phase.residualAlphaRho().value())
            {
                Ei += m[i]/Mvalid*R;
            }
            phase.alphaRhoE()[celli] = Ei;
            phase.alphaRhoU()[celli] = m[i]*VI;
            phase.U()[celli] = VI;
        }
    }

    forAll(fluidPhaseModels_, i)
    {
        phaseModel& phase = fluidPhaseModels_[i];
        phase.U().correctBoundaryConditions();
        phase.alphaRhoU().boundaryFieldRef() ==
            phase.alphaRho().boundaryField()*phase.U().boundaryField();
    }
}


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
        relaxFluidVelocitiesInstant();
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

        volScalarField Kd(dragModelIter()->stabK());

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

            const scalar rAlpha1 = phase1.residualAlpha().value();
            const scalar rAlpha2 = phase2.residualAlpha().value();
            forAll(deltaM, celli)
            {
                if (phase1[celli] <= rAlpha1 || phase2[celli] <= rAlpha2)
                {
                    deltaM[celli] = Zero;
                }
            }

            // Work of the drag force at the interfacial velocity averaged
            // over the implicit update, so that the work done on a phase
            // moving with the interfacial velocity equals the change of its
            // kinetic energy (the dissipation goes to the other phase)
            volVectorField Ui(*Uis[pair]);
            {
                const volVectorField U1(phase1.U());
                const volVectorField U2(phase2.U());
                phase1.U() = U1 + deltaM/alphaRho1;
                phase2.U() = U2 - deltaM/alphaRho2;
                Ui = 0.5*(Ui + interfacialVelocityModels_[pair]->UI());
                phase1.U() = U1;
                phase2.U() = U2;
            }

            phase1.alphaRhoU() += deltaM;
            if (phase1.totalEnergy())
            {
                phase1.alphaRhoE() += deltaM & Ui;
            }

            phase2.alphaRhoU() -= deltaM;
            if (phase2.totalEnergy())
            {
                phase2.alphaRhoE() -= deltaM & Ui;
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

            const scalar rAlphaP = particles->residualAlpha().value();
            const scalar rAlphaG = gas->residualAlpha().value();
            volScalarField deltaPTE
            (
                1.5*particles->alphaRho()*(ThetaStarStar - ThetaOld)
            );
            forAll(deltaPTE, celli)
            {
                if ((*particles)[celli] <= rAlphaP || (*gas)[celli] <= rAlphaG)
                {
                    deltaPTE[celli] = 0.0;
                }
            }

            gas->alphaRhoE() -= deltaPTE;
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

        const scalar rAlpha1 = phase1.residualAlpha().value();
        const scalar rAlpha2 = phase2.residualAlpha().value();
        forAll(deltaE, celli)
        {
            if (phase1[celli] <= rAlpha1 || phase2[celli] <= rAlpha2)
            {
                deltaE[celli] = 0.0;
            }
        }

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
    }
    const volScalarField rhoLim(max(rho_, dimensionedScalar(dimDensity, small)));
    U_ = alphaRhoU/rhoLim;

    forAll(T_, celli)
    {
        scalar activeRho = 0.0;
        scalar activeRhoT = 0.0;
        forAll(phaseModels_, phasei)
        {
            const phaseModel& phase = phaseModels_[phasei];
            const scalar rAlpha = phase.residualAlpha().value();
            if (phase[celli] > rAlpha)
            {
                const scalar aRho = phase.alphaRho()[celli];
                activeRho += aRho;
                activeRhoT += aRho*phase.T()[celli];
            }
        }
        T_[celli] = activeRho > small ? (activeRhoT / activeRho) : 288.15;
    }

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
    sumAlpha.max(1e-12);
    p_ /= sumAlpha;
    sumAlphaRho.max(dimensionedScalar("0", dimDensity, 1e-12));
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
    PRelaxation_(NONE),
    alphaDilute_(lookupOrDefault<scalar>("alphaDilute", 0.01))
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
    //- 1. Decode slave (granular) phases first to calculate compressed density and alpha
    forAll(slavePhaseModels_, i)
    {
        slavePhaseModels_[i].correctBoundaryConditions();
        slavePhaseModels_[i].decode();
    }

    //- 2. Adjust carrier fluid volume fractions so sum(alpha) = 1
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
                const scalar scale = 1.0/fixedAlpha;
                forAll(slavePhaseModels_, i)
                {
                    slavePhaseModels_[i].scaleVolumeFraction
                    (
                        scale,
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

    //- 3. Decode fluid phases now that carrier volume fractions are updated
    forAll(fluidPhaseModels_, i)
    {
        fluidPhaseModels_[i].correctBoundaryConditions();
        fluidPhaseModels_[i].decode();
    }

    //- 4. Pressure equilibrium of dilute fluid phases, then decode again
    if (equilibrateDiluteVolumes())
    {
        forAll(fluidPhaseModels_, i)
        {
            fluidPhaseModels_[i].correctBoundaryConditions();
            fluidPhaseModels_[i].decode();
        }
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


bool Foam::phaseSystem::equilibrateDiluteVolumes()
{
    if (alphaDilute_ <= 0 || fluidPhaseModels_.size() < 2)
    {
        return false;
    }

    UPtrList<fluidBlastThermo> thermos(fluidPhaseModels_.size());
    forAll(fluidPhaseModels_, i)
    {
        thermos.set
        (
            i,
            &mesh_.lookupObjectRef<fluidBlastThermo>
            (
                IOobject::groupName
                (
                    physicalProperties::typeName,
                    fluidPhaseModels_[i].name()
                )
            )
        );
    }

    // Pressure of fluid phase i in celli at the given density and energy
    auto cellp = [&](const label i, const scalar rho, const scalar e, const label celli)
    {
        fluidBlastThermo& thermo = thermos[i];
        thermo.rhoRef()[celli] = rho;
        thermo.he()[celli] = e;
        thermo.T()[celli] = thermo.cellThe(e, thermo.T()[celli], celli);
        return thermo.cellpRhoT(celli);
    };

    // rho*c^2/alpha of fluid phase i at the state last set by cellp
    auto cellC = [&](const label i, const scalar p, const scalar alpha, const label celli)
    {
        fluidBlastThermo& thermo = thermos[i];
        return
            thermo.rhoRef()[celli]
           *sqr(thermo.cellSpeedOfSound(p, celli))
           /alpha;
    };

    const label maxIter = 50;
    const label maxHalvings = 10;
    const scalar pRelTol = 1e-6;
    bool changed = false;

    forAll(rho_, celli)
    {
        // Carrier: the fluid phase with the largest volume fraction
        label c = 0;
        forAll(fluidPhaseModels_, i)
        {
            if (fluidPhaseModels_[i][celli] > fluidPhaseModels_[c][celli])
            {
                c = i;
            }
        }
        phaseModel& carrier = fluidPhaseModels_[c];

        forAll(fluidPhaseModels_, k)
        {
            phaseModel& phase = fluidPhaseModels_[k];
            const scalar alphak0 = phase[celli];
            const scalar alphaRhok = phase.alphaRho()[celli];
            if
            (
                k == c
             || alphak0 >= alphaDilute_
             || alphak0 <= phase.residualAlpha().value()
             || alphaRhok <= phase.residualAlphaRho().value()
            )
            {
                continue;
            }

            const scalar alphac0 = carrier[celli];
            const scalar alphaRhoc = carrier.alphaRho()[celli];
            const scalar Tk0 = thermos[k].T()[celli];
            const scalar Tc0 = thermos[c].T()[celli];
            // Energy bounds of the thermo temperature range
            const scalar ekMin =
                thermos[k].cellhe(thermos[k].TLow(), celli);
            const scalar ekMax =
                thermos[k].cellhe(thermos[k].THigh(), celli);
            const scalar ecMin =
                thermos[c].cellhe(thermos[c].TLow(), celli);
            const scalar ecMax =
                thermos[c].cellhe(thermos[c].THigh(), celli);

            // Damped Newton iteration on the volume exchanged; the work of
            // each increment is evaluated at the mean pressure, so total
            // energy is conserved
            scalar alphak = alphak0;
            scalar alphac = alphac0;
            const scalar ek0 = thermos[k].he()[celli];
            const scalar ec0 = thermos[c].he()[celli];
            scalar ek = ek0;
            scalar ec = ec0;
            scalar pk = cellp(k, alphaRhok/alphak, ek, celli);
            scalar Ck = cellC(k, pk, alphak, celli);
            scalar pc = cellp(c, alphaRhoc/alphac, ec, celli);
            scalar Cc = cellC(c, pc, alphac, celli);

            for (label iter = 0; iter < maxIter; iter++)
            {
                const scalar dp = pk - pc;
                if (mag(dp) <= pRelTol*(mag(pk) + mag(pc)) + small)
                {
                    break;
                }

                scalar dAlpha = dp/max(Ck + Cc, small);

                // At most double the dilute phase volume and take at most
                // 10% of the carrier volume per iteration. The dilute phase
                // is compressed by at most a factor of two per call:
                // compressing e.g. trace gas in a detonating explosive to
                // the carrier pressure can take it far outside the range of
                // its equation of state.
                dAlpha = min(dAlpha, alphak);
                dAlpha = min(dAlpha, 0.1*alphac);
                dAlpha = max(dAlpha, 0.5*alphak0 - alphak);
                if (dAlpha == 0)
                {
                    break;
                }

                const scalar pI = 0.5*(pk + pc);
                bool accepted = false;
                scalar ekNew = ek, ecNew = ec, pkNew = pk, pcNew = pc;
                for (label h = 0; h < maxHalvings; h++)
                {
                    ekNew = ek - pI*dAlpha/alphaRhok;
                    ecNew = ec + pI*dAlpha/alphaRhoc;
                    if
                    (
                        ekNew > ekMin && ekNew < ekMax
                     && ecNew > ecMin && ecNew < ecMax
                    )
                    {
                        pkNew = cellp(k, alphaRhok/(alphak + dAlpha), ekNew, celli);
                        pcNew = cellp(c, alphaRhoc/(alphac - dAlpha), ecNew, celli);
                        if (mag(pkNew - pcNew) < mag(dp))
                        {
                            accepted = true;
                            break;
                        }
                    }
                    dAlpha *= 0.5;
                }
                if (!accepted)
                {
                    break;
                }

                alphak += dAlpha;
                alphac -= dAlpha;
                ek = ekNew;
                ec = ecNew;
                pk = pkNew;
                pc = pcNew;
                Ck = cellC(k, pk, alphak, celli);
                Cc = cellC(c, pc, alphac, celli);
            }

            // Restore the thermo temperatures (the cell is decoded again)
            thermos[k].T()[celli] = Tk0;
            thermos[c].T()[celli] = Tc0;

            if (alphak == alphak0)
            {
                continue;
            }

            phase[celli] = alphak;
            carrier[celli] = alphac;
            phase.alphaRhoE()[celli] += alphaRhok*(ek - ek0);
            carrier.alphaRhoE()[celli] += alphaRhoc*(ec - ec0);
            changed = true;
        }
    }

    return returnReduce(changed, orOp<bool>());
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
        relaxFluidVelocitiesInstant();
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


void Foam::phaseSystem::solveExplicit()
{
    decode();

    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].solveExplicit();
    }

    const dimensionedScalar& deltaT(mesh_.time().deltaT());
    relaxVelocity(deltaT);
    relaxTemperature(deltaT);

    master_.solveExplicit();
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


void Foam::phaseSystem::solveImplicit()
{
    decode();

    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].solveImplicit();
    }
    master_.solveImplicit();

    relaxPressure(mesh_.time().deltaT());
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


void Foam::phaseSystem::postUpdate()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].postUpdate();
    }
    master_.postUpdate();
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
        mD21*(phase2.thermo().hs()()() + (hc2()() - hc1()()))
      + mD12*phase1.thermo().hs()()();

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
