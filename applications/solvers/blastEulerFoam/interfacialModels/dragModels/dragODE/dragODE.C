/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
7-2-2019 Jeff Heylmun:      Added drag terms
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

#include "dragODE.H"

Foam::label Foam::dragODE::calcNEqns()
{
    label i = 0;
    forAll(phaseModels_, phasei)
    {
        startI_[phasei] = i;
        i += nDims_*phaseModels_[phasei].nNodes();
    }
    return i;
}

void Foam::dragODE::setUs()
{
    forAll(phaseModels_, phasei)
    {
        phaseModel& phase = phaseModels_[phasei];
        for (label nodei = 0; nodei < phase.nNodes(); nodei ++)
        {
            label ui = startI_[phasei];
            if (solutionD_[0] == 1)
            {
                phase.U(nodei)[celli_].x() = q_[ui];
                ui++;
            }
            if (solutionD_[1] == 1)
            {
                phase.U(nodei)[celli_].y() = q_[ui];
                ui++;

            }
            if (solutionD_[2] == 1)
            {
                phase.U(nodei)[celli_].z() = q_[ui];
            }
        }
    }
}


void Foam::dragODE::setq()
{
    forAll(phaseModels_, phasei)
    {
        phaseModel& phase = phaseModels_[phasei];
        label ui = startI_[phasei];
        for (label nodei = 0; nodei < phase.nNodes(); nodei++)
        {
            vector U = phase.U(nodei)[celli_];
            if (solutionD_[0] == 1)
            {
                q_[ui] = U.x();
                ui++;
            }
            if (solutionD_[1] == 1)
            {
                q_[ui] = U.y();
                ui++;
            }
            if (solutionD_[2] == 1)
            {
                q_[ui] = U.z();
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragODE::dragODE(phaseSystem& fluid, dragModelTable& dragModels)
:
    ODESystem(),
    dict_(fluid.subDict("dragODECoeffs")),
    solveDrag_(dict_.lookupOrDefault("solveODE", false)),
    fluid_(fluid),
    phaseModels_(fluid.phases()),
    dragModels_(dragModels.size()),
    phasePairs_(dragModels.size()),
    solutionD_(fluid.mesh().solutionD()),
    nDims_(fluid_.mesh().nSolutionD()),
    startI_(phaseModels_.size()),
    nEqns_(calcNEqns()),
    q_(nEqns_, 0.0),
    dqdt_(nEqns_, 0.0),
    deltaT_(phaseModels_[0].size(), fluid.mesh().time().deltaTValue())
{
    odeSolver_ = ODESolver::New(*this, dict_);

    //- Add unorded phase pairs with vaild drag models
    label i = 0;
    forAllIter
    (
        dragModelTable,
        dragModels,
        dragModelIter
    )
    {
        const phasePair& pair(fluid.phasePairs()[dragModelIter.key()]);
        const phasePairKey key(pair.first(), pair.second());

        phasePairs_.set(i, new phasePair(pair));
        dragModels_.set(i, &dragModelIter()());
        i++;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragODE::~dragODE()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::dragODE::derivatives
(
    const scalar time,
    const scalarField& q,
    scalarField& dqdt
) const
{
    dqdt = 0.0;

    forAll(phasePairs_, pairi)
    {
        const phasePair& pair(phasePairs_[pairi]);
        const phaseModel& phase1(pair.phase1());
        const phaseModel& phase2(pair.phase2());
        const scalar rAlphaRho1 = phase1.residualAlphaRho().value();
        const scalar rAlphaRho2 = phase2.residualAlphaRho().value();
        const scalar rho1 = phase1.rho()[celli_];
        const scalar rho2 = phase2.rho()[celli_];

        for (label nodei = 0; nodei < phase1.nNodes(); nodei++)
        {
            label ui = startI_[phase1.index()];
            scalar alphaRho1 =
                Foam::max
                (
                    phase1.volumeFractioni(celli_, nodei)*rho1,
                    rAlphaRho1
                );
            for (label nodej = 0; nodej < phase2.nNodes(); nodej++)
            {

                label uj = startI_[phase2.index()];
                scalar alphaRho2 =
                    Foam::max
                    (
                        phase2.volumeFractioni(celli_, nodei)*rho2,
                        rAlphaRho2
                    );

                scalar drag = dragModels_[pairi].Ki(celli_, nodei, nodej);
                scalar drag1 = drag/alphaRho1;
                scalar drag2 = drag/alphaRho2;

                if (nDims_ > 0)
                {
                    scalar u1 = q[ui];
                    scalar u2 = q[uj];

                    dqdt[ui] += drag1*(u2 - u1);
                    dqdt[uj] += drag2*(u1 - u2);
                    ui++;
                    uj++;
                }
                if (nDims_ > 1)
                {
                    scalar u1 = q[ui];
                    scalar u2 = q[uj];

                    dqdt[ui] += drag1*(u2 - u1);
                    dqdt[uj] += drag2*(u1 - u2);
                    ui++;
                    uj++;
                }
                if (nDims_ > 2)
                {
                    scalar u1 = q[ui];
                    scalar u2 = q[uj];
                    dqdt[ui] += drag1*(u2 - u1);
                    dqdt[uj] += drag2*(u1 - u2);
                }
            }
        }
    }

}


void Foam::dragODE::jacobian
(
    const scalar t,
    const scalarField& q,
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
        const scalar rAlpha1 = phase1.residualAlpha().value();
        const scalar rAlpha2 = phase2.residualAlpha().value();
        const scalar rho1 = phase1.rho()[celli_];
        const scalar rho2 = phase2.rho()[celli_];

        for (label nodei = 0; nodei < phase1.nNodes(); nodei++)
        {
            label ui = startI_[phase1.index()];
            scalar alphaRho1 =
                Foam::max
                (
                    phase1.volumeFractioni(celli_, nodei),
                    rAlpha1
                )*rho1;
            for (label nodej = 0; nodej < phase2.nNodes(); nodej++)
            {

                label uj = startI_[phase2.index()];
                scalar alphaRho2 =
                    Foam::max
                    (
                        phase2.volumeFractioni(celli_, nodei),
                        rAlpha2
                    )*rho2;

                scalar drag = dragModels_[pairi].Ki(celli_, nodei, nodej);
                scalar drag1 = drag/alphaRho1;
                scalar drag2 = drag/alphaRho2;

                if (nDims_ > 0)
                {
                    scalar u1 = q[ui];
                    scalar u2 = q[uj];

                    dqdt[ui] += drag1*(u2 - u1);
                    dqdt[uj] += drag2*(u1 - u2);

                    J(ui, ui) -= drag1;
                    J(ui, uj) += drag1;
                    J(uj, uj) -= drag2;
                    J(uj, ui) += drag2;

                    ui++;
                    uj++;
                }
                if (nDims_ > 1)
                {
                    scalar u1 = q[ui];
                    scalar u2 = q[uj];

                    dqdt[ui] += drag1*(u2 - u1);
                    dqdt[uj] += drag2*(u1 - u2);

                    J(ui, ui) -= drag1;
                    J(ui, uj) += drag1;
                    J(uj, uj) -= drag2;
                    J(uj, ui) += drag2;

                    ui++;
                    uj++;
                }
                if (nDims_ > 2)
                {
                    scalar u1 = q[ui];
                    scalar u2 = q[uj];

                    dqdt[ui] += drag1*(u2 - u1);
                    dqdt[uj] += drag2*(u1 - u2);
                    J(ui, ui) -= drag1;
                    J(ui, uj) += drag1;
                    J(uj, uj) -= drag2;
                    J(uj, ui) += drag2;
                }
            }
        }
    }
}


Foam::scalar Foam::dragODE::solve
(
    const scalar& deltaT
)
{
    for(celli_ = 0; celli_ < phaseModels_[0].size(); celli_++)
    {
        scalarList oldKs(phaseModels_.size());
        forAll(oldKs, phasei)
        {
            oldKs[phasei] =
                phaseModels_[phasei].alphaRho()[celli_]
               *0.5*magSqr(phaseModels_[phasei].U()[celli_]);
        }
        setq();

        scalar timeLeft = deltaT;
        while (timeLeft > small)
        {
            scalar dt = timeLeft;
            odeSolver_->solve(0, dt, q_, deltaT_[celli_]);
            setUs();
            timeLeft -= dt;
            deltaT_[celli_] = dt;
        }

        forAll(phaseModels_, phasei)
        {
            phaseModels_[phasei].alphaRhoE()[celli_] +=
                phaseModels_[phasei].alphaRho()[celli_]
               *0.5*magSqr(phaseModels_[phasei].U()[celli_])
              - oldKs[phasei];
        }
    }
    return min(deltaT_);
}


// ************************************************************************* //
