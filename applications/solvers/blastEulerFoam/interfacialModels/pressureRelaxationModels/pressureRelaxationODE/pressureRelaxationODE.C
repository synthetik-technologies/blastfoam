/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pressureRelaxationODE, 0);
    addToRunTimeSelectionTable
    (
        pressureRelaxationSolver,
        pressureRelaxationODE,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressureRelaxationODE::pressureRelaxationODE
(
    phaseSystem& fluid,
    interfacialPressureModelTable& interfacialPressureModels,
    pressureRelaxationModelTable& pressureRelaxationModels
)
:
    pressureRelaxationSolver
    (
        fluid,
        interfacialPressureModels,
        pressureRelaxationModels,
        true,
        true
    ),
    ODESystem(),
    dict_(fluid.subDict("pressureSolverCoeffs")),
    q_(),
    dqdt_(),
    deltaT_
    (
        IOobject
        (
            typeName + ":deltaT",
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        fluid.mesh().time().deltaT()
    )
{
    nEqns_ = phaseModels_.size()*2;
    q_.resize(nEqns_);
    dqdt_.resize(nEqns_);

    if (pressureRelaxationModels_.size() != interfacialPressureModels_.size())
    {
        FatalErrorInFunction
            << "Different number of pressureRelaxationModels and "
            << "interfacialPressureModels" << endl
            << abort(FatalError);
    }
    forAll(pressureRelaxationModels_, i)
    {
        if
        (
            &pressureRelaxationModels_[i].pair()
         != &interfacialPressureModels_[i].pair()
        )
        {

            FatalErrorInFunction
                << "Mismatched pairs" << endl
                << abort(FatalError);
        }
    }

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
    forAll(thermos_, phasei)
    {
        thermos_[phasei].rho()[li] =
            phaseModels_[phasei].alphaRho()[li]
           /max(phaseModels_[phasei][li], 1e-10);
    }
    forAll(interfacialPressureModels_, pairi)
    {
        const interfacialPressureModel& ip =
            interfacialPressureModels_[pairi];
        const phasePair& pair(ip.pair());
        const phaseModel& phase1(pair.phase1());
        const phaseModel& phase2(pair.phase2());

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
        const scalar mu = pressureRelaxationModels_[pairi].cellK(li);

        dqdt[a1i] += mu*(p1 - p2);
        dqdt[e1i] += PI*mu*(p2 - p1)/phase1.alphaRho()[li];
        dqdt[a2i] += mu*(p2 - p1);
        dqdt[e2i] += PI*mu*(p1 - p2)/phase2.alphaRho()[li];
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
    forAll(thermos_, phasei)
    {
        thermos_[phasei].rho()[li] =
            phaseModels_[phasei].alphaRho()[li]
           /max(phaseModels_[phasei][li], 1e-10);
    }
    forAll(interfacialPressureModels_, pairi)
    {
        const interfacialPressureModel& ip =
            interfacialPressureModels_[pairi];
        const phasePair& pair(ip.pair());
        const phaseModel& phase1(pair.phase1());
        const phaseModel& phase2(pair.phase2());

        const label a1i = phaseIndicies_[phase1.index()];
        const label e1i = phaseModels_.size() + a1i;
        const label a2i = phaseIndicies_[phase2.index()];
        const label e2i = phaseModels_.size() + a2i;
        const scalar p1 = thermos_[a1i].cellpRhoT(li);
        const scalar p2 = thermos_[a2i].cellpRhoT(li);

        if (phase1[li] < 1e-10 || phase2[li] < 1e-10)
        {
            continue;
        }

        const scalar alphaRho1 = phase1.alphaRho()[li];
        const scalar alphaRho2 = phase2.alphaRho()[li];

        const scalar dPIdAlpha1 = ip.celldPIdAlpha(li, phase1.index());
        const scalar dPIdAlpha2 = ip.celldPIdAlpha(li, phase2.index());
        const scalar dPIde1 = ip.celldPIde(li, phase1.index());
        const scalar dPIde2 = ip.celldPIde(li, phase2.index());
        const scalar dp1de1 = thermos_[a1i].celldpde(li);
        const scalar dp2de2 = thermos_[a2i].celldpde(li);

        const scalar PI = ip.cellPI(li);
        const scalar mu = pressureRelaxationModels_[pairi].cellK(li);

        dqdt[a1i] += mu*(p1 - p2);
        dqdt[e1i] += PI*mu*(p2 - p1)/alphaRho1;
        dqdt[a2i] += mu*(p2 - p1);
        dqdt[e2i] += PI*mu*(p1 - p2)/alphaRho2;

//         J[a1i][a1i] += 0.0;
        J[e1i][a1i] += mu*dPIdAlpha1*(p2 - p1)/alphaRho1;

//         J[a1i][a2i] += 0.0;
        J[e1i][a2i] += mu*dPIdAlpha2*(p2 - p1)/alphaRho1;

        J[a1i][e1i] += mu*dp1de1;
        J[e1i][e1i] += mu*(dPIde1*(p2 - p1) - PI*dp1de1)/alphaRho1;

        J[a1i][e2i] -= mu*dp2de2;
        J[e1i][e2i] += mu*(dPIde2*(p2 - p1) + PI*dp2de2)/alphaRho1;

//         J[a2i][a2i] -= 0.0;
        J[e2i][a2i] += mu*dPIdAlpha2*(p1 - p2)/alphaRho2;

//         J[a2i][a1i] -= 0.0;
        J[e2i][a1i] += dPIdAlpha1*mu*(p1 - p2)/alphaRho2;

        J[a2i][e2i] += mu*dp2de2;
        J[e2i][e2i] += mu*(dPIde2*(p1 - p2) - PI*dp2de2)/alphaRho2;

        J[a2i][e1i] -= mu*dp1de1;
        J[e2i][e1i] += mu*(dPIde1*(p1 - p2) + PI*dp1de1)/alphaRho2;
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
            const label ai = phaseIndicies_[phasei];
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
    forAll(thermos_, phasei)
    {
        thermos_[phasei].correct();
    }
    return min(deltaT_).value();
}


// ************************************************************************* //
