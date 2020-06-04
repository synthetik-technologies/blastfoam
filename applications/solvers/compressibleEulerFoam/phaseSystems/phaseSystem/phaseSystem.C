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
#include "dragODE.H"
#include "surfaceInterpolate.H"
#include "fvcDdt.H"

#include "dragModel.H"

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
    if (dragOde_.valid())
    {
        Info<< "Solving drag ODE system" <<endl;
        dragOde_->solve(deltaT.value());
        forAll(phaseModels_, phasei)
        {
            phaseModel& phase(phaseModels_[phasei]);
            phase.alphaRhoU() = phase.alphaRho()*phase.U();
        }
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

                if (!dragOde_.valid())
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

                    phase1.alphaRhoU(nodei) += deltaM;
                    if (!phase1.granular())
                    {
                        phase1.alphaRhoE() += deltaM & U_;
                    }
                    phase2.alphaRhoU(nodej) -= deltaM;
                    if (!phase2.granular())
                    {
                        phase2.alphaRhoE() -= deltaM & U_;
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

                volVectorField liftForce
                (
                    liftModels_[pair]->F(nodei, nodej)*deltaT
                );

                phase1.alphaRhoU(nodei) += liftForce;
                if (!phase1.granular())
                {
                    phase1.alphaRhoE() += liftForce & U_;
                }
                phase2.alphaRhoU(nodej) -= liftForce;
                if (!phase2.granular())
                {
                    phase2.alphaRhoE() -= liftForce & U_;
                }
            }
        }
    }

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


void Foam::phaseSystem::calcRho()
{
    const label nPhases = phaseModels_.size();

    rho_ = phaseModels_[0].alphaRho();
    for (label phasei = 1; phasei < nPhases; ++ phasei)
    {
        rho_ += phaseModels_[phasei].alphaRho();
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

    phaseModels_(lookup("phases"), phaseModel::iNew(*this)),

    interfacialPressure_
    (
        interfacialPressureModel::New
        (
            subDict("interfacialPressure"),
            phaseModels_
        )
    ),
    interfacialVelocity_
    (
        interfacialVelocityModel::New
        (
            subDict("interfacialVelocity"),
            phaseModels_
        )
    ),

    kineticTheoryPtr_(NULL),
    polydisperseKineticTheory_(false)
{
    // Sub-models
    generatePairsAndSubModels("aspectRatio", aspectRatioModels_);
    generatePairsAndSubModels("drag", dragModels_);
    generatePairsAndSubModels("virtualMass", virtualMassModels_);
    generatePairsAndSubModels("lift", liftModels_);
    generatePairsAndSubModels("turbulentDispersion", turbulentDispersionModels_);
    generatePairsAndSubModels("wallLubrication", wallLubricationModels_);
    generatePairsAndSubModels("heatTransfer", heatTransferModels_);

    if (lookupOrDefault<Switch>("solveDragODE", false))
    {
        dragOde_.set(new dragOde(*this, dragModels_));
    }

    phaseModels_[phaseModels_.size() - 1].solveAlpha(false);
    volScalarField sumAlpha(phaseModels_[0]);
    label nPhases = phaseModels_.size();
    for (label phasei = 1; phasei < nPhases - 1; phasei++)
    {
        sumAlpha += phaseModels_[phasei];
    }
    sumAlpha.min(1.0);
    volScalarField& alpha(phaseModels_[nPhases - 1]);
    alpha = 1.0 - sumAlpha;
    alpha.correctBoundaryConditions();
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

    calcRho();
    U_ = interfacialVelocity_->Ui();
    p_ = interfacialPressure_->Pi();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseSystem::~phaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseSystem::decode()
{
    phaseModels_[0].decode();
    volScalarField sumAlpha(phaseModels_[0]);
    label nPhases = phaseModels_.size();
    for (label phasei = 1; phasei < nPhases - 1; phasei++)
    {
        phaseModels_[phasei].decode();
        sumAlpha += phaseModels_[phasei];
    }
    sumAlpha.min(1.0);
    volScalarField& alpha(phaseModels_[nPhases - 1]);
    alpha = 1.0 - sumAlpha;
    alpha.max(0);
    alpha.min(1.0);
    alpha.correctBoundaryConditions();
    phaseModels_[nPhases - 1].decode();


    if (kineticTheoryPtr_ != NULL)
    {
        kineticTheoryPtr_->correct();
    }

    calcRho();
    U_ = interfacialVelocity_->Ui();
    p_ = interfacialPressure_->Pi();
}


Foam::tmp<Foam::surfaceScalarField> Foam::phaseSystem::phi() const
{
    return interfacialVelocity_->phi();
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

    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].update();
    }
}


void Foam::phaseSystem::solve
(
    const label stepi,
    const scalarList& ai,
    const scalarList& bi
)
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].solve(stepi, ai, bi);
    }
    decode();
}


void Foam::phaseSystem::setODEFields
(
    const label nSteps,
    const boolList& storeFields,
    const boolList& storeDeltas
)
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].setODEFields(nSteps, storeFields, storeDeltas);
    }
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


void Foam::phaseSystem::relax()
{
    const dimensionedScalar& deltaT(mesh_.time().deltaT());

    relaxVelocity(deltaT);
    relaxTemperature(deltaT);

    decode();
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
