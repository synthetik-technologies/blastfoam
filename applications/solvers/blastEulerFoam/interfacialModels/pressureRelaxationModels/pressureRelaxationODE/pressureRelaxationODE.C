/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
7-2-2019 Jeff Heylmun:      Added pressureRelaxation terms
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

#include "pressureRelaxationODE.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressureRelaxationODE::pressureRelaxationODE(phaseSystem& fluid, pressureRelaxationModelTable& pressureRelaxationModels)
:
    ODESystem(),
    dict_(fluid.subDict("pressureODECoeffs")),
    solvePressureRelaxation_(dict_.lookupOrDefault("solveODE", true)),
    fluid_(fluid),
    phaseModels_(0),
    phaseIndicies_(fluid.phases().size(), -1),
    thermos_(0),
    interfacialPressureModels_(pressureRelaxationModels.size()),
    pressureRelaxationModels_(pressureRelaxationModels.size()),
    phasePairs_(pressureRelaxationModels.size()),
    nEqns_(0),
    q_(0),
    dqdt_(0),
    deltaT_
    (
        IOobject
        (
            "pressureRelaxationODE",
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        fluid.mesh().time().deltaT()
    )
{
    label nFluids = 0;
    forAll(fluid.phases(), phasei)
    {
        if (!fluid.phases()[phasei].slavePressure())
        {
            phaseIndicies_[phasei] = nFluids;
            nFluids++;
        }
    }
    nEqns_ = 2*nFluids; // alpha and internal energy
    q_.resize(nEqns_);
    dqdt_.resize(nEqns_);

    if (nFluids <= 1)
    {
        solvePressureRelaxation_ = false;
    }

    phaseModels_.resize(nFluids);
    thermos_.resize(nFluids);
    label i = 0;
    forAll(fluid.phases(), phasei)
    {
        const phaseModel& phase = fluid.phases()[phasei];
        if (!phase.slavePressure())
        {
            phaseModels_.set
            (
                i,
                &fluid.mesh().lookupObjectRef<phaseModel>
                (
                    IOobject::groupName("alpha", phase.name())
                )
            );
            thermos_.set
            (
                i,
                &fluid.mesh().lookupObjectRef<fluidBlastThermo>
                (
                    IOobject::groupName(basicThermo::dictName, phase.group())
                )
            );
            i++;
        }
    }

    //- Add unorded phase pairs with vaild pressureRelaxation models
    i = 0;

    forAllIter
    (
        pressureRelaxationModelTable,
        pressureRelaxationModels,
        pressureRelaxationModelIter
    )
    {
        const phasePair& pair(fluid.phasePairs()[pressureRelaxationModelIter.key()]);
        const phasePairKey key(pair.first(), pair.second());

        phasePairs_.set(i, new phasePair(pair));
        interfacialPressureModels_.set
        (
            i,
            &const_cast<interfacialPressureModel&>
            (
                fluid.lookupSubModel<interfacialPressureModel>(pair)
            )
        );
        pressureRelaxationModels_.set(i, &pressureRelaxationModelIter()());
        i++;
    }
    phasePairs_.resize(i);
    interfacialPressureModels_.resize(i);
    pressureRelaxationModels_.resize(i);

    odeSolver_ = ODESolver::New(*this, dict_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pressureRelaxationODE::~pressureRelaxationODE()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::pressureRelaxationODE::derivatives
(
    const scalar time,
    const scalarField& q,
    const label li,
    scalarField& dqdt
) const
{
    dqdt = scalarField(nEqns_, 0.0);
    forAll(phasePairs_, pairi)
    {
        const phasePair& pair(phasePairs_[pairi]);
        const phaseModel& phase1(pair.phase1());
        const phaseModel& phase2(pair.phase2());
        const interfacialPressureModel& ip =
            interfacialPressureModels_[pairi];

        const label a1i = phaseIndicies_[phase1.index()];
        const label e1i = phaseModels_.size() + a1i;
        const label a2i = phaseIndicies_[phase2.index()];
        const label e2i = phaseModels_.size() + a2i;
        const scalar p1 = thermos_[a1i].cellpRhoT(li);
        const scalar p2 = thermos_[a2i].cellpRhoT(li);
        if
        (
            phase1.alphaRho()[li] < 1e-10
         || phase2.alphaRho()[li] < 1e-10
        )
        {
            continue;
        }
        const scalar PI = ip.cellPI(li);

        const scalar alphaRho1 = phase1.alphaRho()[li];
        const scalar alphaRho2 = phase2.alphaRho()[li];

        const scalar mu = pressureRelaxationModels_[pairi].cellK(li);

        dqdt[a1i] += mu*(p1 - p2);
        dqdt[e1i] += PI*mu*(p2 - p1)/alphaRho1;
        dqdt[a2i] -= mu*(p1 - p2);
        dqdt[e2i] -= PI*mu*(p2 - p1)/alphaRho2;
    }
}


void Foam::pressureRelaxationODE::jacobian
(
    const scalar t,
    const scalarField& q,
    const label li,
    scalarField& dqdt,
    scalarSquareMatrix& J
) const
{
    dqdt = 0.0;
    J = scalarSquareMatrix(nEqns_, 0.0);
    forAll(phasePairs_, pairi)
    {
        const phasePair& pair(phasePairs_[pairi]);
        const phaseModel& phase1(pair.phase1());
        const phaseModel& phase2(pair.phase2());
        const interfacialPressureModel& ip =
            interfacialPressureModels_[pairi];
        const label a1i = phaseIndicies_[phase1.index()];
        const label e1i = phaseModels_.size() + a1i;
        const label a2i = phaseIndicies_[phase2.index()];
        const label e2i = phaseModels_.size() + a2i;
        const scalar p1 = thermos_[a1i].cellpRhoT(li);
        const scalar p2 = thermos_[a2i].cellpRhoT(li);

        const scalar& rho1 = phase1.rho()[li];
        const scalar& rho2 = phase2.rho()[li];

        if
        (
            phase1.alphaRho()[li] < 1e-10
         || phase2.alphaRho()[li] < 1e-10
        )
        {
            continue;
        }
        const scalar PI = ip.cellPI(li);

        const scalar alphaRho1 = phase1.alphaRho()[li];
        const scalar alphaRho2 = phase2.alphaRho()[li];

        const scalar dPidAlpha1 = ip.celldPIdAlpha(li, phase1.index());
        const scalar dPidAlpha2 = ip.celldPIdAlpha(li, phase2.index());
        const scalar dPide1 = ip.celldPIde(li, phase1.index());
        const scalar dPide2 = ip.celldPIde(li, phase2.index());
        const scalar dp1de1 = thermos_[a1i].celldpde(li);
        const scalar dp2de2 = thermos_[a2i].celldpde(li);

        const scalar mu = pressureRelaxationModels_[pairi].cellK(li);

        dqdt[a1i] += mu*(p1 - p2);
        dqdt[e1i] += PI*mu*(p2 - p1)/alphaRho1;
        dqdt[a2i] -= mu*(p1 - p2);
        dqdt[e2i] -= PI*mu*(p2 - p1)/alphaRho2;

        J[a1i][a1i] += 0.0;
        J[e1i][a1i] +=
            (dPidAlpha1 - rho1*PI/alphaRho1)*mu*(p2 - p1)/alphaRho1;

        J[a1i][a2i] += 0.0;
        J[e1i][a2i] += dPidAlpha2*mu*(p2 - p1)/alphaRho1;

        J[a1i][e1i] += mu*dPide1;
        J[e1i][e1i] += mu/alphaRho1*((p2 - p1)*dPide1 - PI*dp1de1);

        J[a1i][e2i] += -mu*dPide2;
        J[e1i][e2i] += mu/alphaRho1*((p2 - p1)*dPide2 + PI*dp2de2);

        J[a2i][a2i] -= 0.0;
        J[e2i][a2i] -=
            (dPidAlpha2 - rho2*PI/alphaRho2)*mu*(p2 - p1)/alphaRho2;

        J[a2i][a1i] -= 0.0;
        J[e2i][a1i] -= dPidAlpha1*mu*(p2 - p1)/alphaRho2;

        J[a2i][e2i] -= mu*dPide2;
        J[e2i][e2i] -= mu/alphaRho2*((p2 - p1)*dPide2 + PI*dp2de2);

        J[a2i][e1i] -= -mu*dPide1;
        J[e2i][e1i] -= mu/alphaRho2*((p2 - p1)*dPide1 - PI*dp1de1);
    }
}


Foam::scalar Foam::pressureRelaxationODE::solve
(
    const scalar& deltaT
)
{
    forAll(phaseModels_[0], celli)
    {
        forAll(phaseModels_, phasei)
        {
            const phaseModel& phase = phaseModels_[phasei];
            const label ai = phaseIndicies_[phase.index()];
            const label ei = phaseModels_.size() + ai;
            q_[ai] = max(min(phase[celli], 1.0), 0.0);
            q_[ei] = phase.he()[celli];
        }

        scalar timeLeft = deltaT;
        while (timeLeft > small)
        {
            scalar dt = timeLeft;
            odeSolver_->solve(0, dt, q_, celli, deltaT_[celli]);
            timeLeft -= dt;
            deltaT_[celli] = dt;

            forAll(phaseModels_, phasei)
            {
                phaseModel& phase = phaseModels_[phasei];
                const label ai = phaseIndicies_[phase.index()];
                const label ei = phaseModels_.size() + ai;
                phase[celli] = max(min(q_[ai], 1.0), 0.0);
                phase.rho()[celli] =
                    phase.alphaRho()[celli]/max(phase[celli], 1e-10);
                phase.he()[celli] = q_[ei];
            }
        }

    }
    return min(deltaT_).value();
}


// ************************************************************************* //
