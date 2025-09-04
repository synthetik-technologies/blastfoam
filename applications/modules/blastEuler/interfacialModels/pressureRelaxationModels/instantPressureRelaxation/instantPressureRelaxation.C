/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022-2025
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

#include "instantPressureRelaxation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(instantPressureRelaxation, 0);
    addToRunTimeSelectionTable
    (
        pressureRelaxationSolver,
        instantPressureRelaxation,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::instantPressureRelaxation::instantPressureRelaxation
(
    const dictionary& dict,
    phaseSystem& fluid,
    interfacialPressureModelTable& interfacialPressureModels,
    pressureRelaxationModelTable& pressureRelaxationModels
)
:
    pressureRelaxationSolver(fluid, true),
    ScalarMultivariateEquation
    (
        phaseModels_.size() + 1,
        scalarField(phaseModels_.size() + 1, 0.0),
        scalarField(phaseModels_.size() + 1, great)
    ),
    useRootSolver_
    (
        phaseModels_.size() > 2
     || dict.lookupOrDefault<bool>("useRootSolver", true)
    )
{
    if (useRootSolver_)
    {
        // Set dX values
        this->dX_ = 1e-6;
        this->dX_.last() = 1.0;

        // Set number of equation (alpha1->alphaN, PI)
        pressureRelaxationSolver::nEqns_ = phaseModels_.size() + 1;

        // No upper limit for any
        rootSolver_ = rootSolver::New(*this, dict);

        rootSolver_->relTolerances() = 1e-10;
        rootSolver_->relTolerances().last() = 1e-10;
        rootSolver_->absTolerances() = small;
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::instantPressureRelaxation::~instantPressureRelaxation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::instantPressureRelaxation::relax4
(
    const scalar deltaT
)
{
    const label maxIter = 100;
    const scalar pRelTol = 1e-6;
    phaseModel& phase1 = phaseModels_[0];
    phaseModel& phase2 = phaseModels_[1];

    fluidBlastThermo& thermo1 = thermos_[0];
    fluidBlastThermo& thermo2 = thermos_[1];

    const scalar rAlpha1 = phase1.residualAlpha().value();
    const scalar rAlpha2 = phase2.residualAlpha().value();

    const scalar rAlphaRho1 = phase1.residualAlphaRho().value();
    const scalar rAlphaRho2 = phase2.residualAlphaRho().value();

    forAll(phase1, celli)
    {
        scalar& alpha1 = phase1[celli];
        scalar& alpha2 = phase2[celli];

        const scalar alphaRho1 = phase1.alphaRho()[celli];
        const scalar alphaRho2 = phase2.alphaRho()[celli];

        if
        (
            alpha1 < rAlpha1
         || alpha2 < rAlpha2
         || alphaRho1 < rAlphaRho1
         || alphaRho2 < rAlphaRho2
        )
        {
            continue;
        }

        scalar& rho1 = thermo1.rhoRef()[celli];
        scalar& rho2 = thermo2.rhoRef()[celli];

        rho1 = alphaRho1/max(alpha1, small);
        rho2 = alphaRho2/max(alpha2, small);

        scalar& e1 = thermo1.he()[celli];
        scalar& e2 = thermo2.he()[celli];

        const scalar KE1 = 0.5*magSqr(phase1.alphaRhoU()[celli]/alphaRho1);
        const scalar KE2 = 0.5*magSqr(phase2.alphaRhoU()[celli]/alphaRho2);

        e1 = phase1.alphaRhoE()[celli]/alphaRho1 - KE1;
        e2 = phase2.alphaRhoE()[celli]/alphaRho2 - KE2;

        scalar& T1 = thermo1.T()[celli];
        scalar& T2 = thermo2.T()[celli];

        scalar e1Test = thermo1.cellhe(thermo1.TLow(), celli);
        if (e1 < e1Test)
        {
            e1 = e1Test;
            T1 = thermo1.TLow();
        }
        else
        {
            T1 = thermo1.cellThe(e1, T1, celli);
        }

        scalar e2Test = thermo2.cellhe(thermo2.TLow(), celli);
        if (e2 < e2Test)
        {
            e2 = e2Test;
            T2 = thermo2.TLow();
        }
        else
        {
            T2 = thermo2.cellThe(e2, T2, celli);
        }

        scalar& p1 = thermo1.p()[celli];
        scalar& p2 = thermo2.p()[celli];

        p1 = phase1.cellpRhoT(celli);
        p2 = phase2.cellpRhoT(celli);

        scalar PI = fluid_.cellPI(celli);

        scalar C1 = sqr(phase1.cellSpeedOfSound(PI, celli))*rho1/alpha1;
        scalar C2 = sqr(phase2.cellSpeedOfSound(PI, celli))*rho2/alpha2;

        scalar dPI = p1 - p2;
        scalar dPIOld = dPI;

        if (mag(dPI)/PI < pRelTol)
        {
            continue;
        }

        label iter = 0;
        do
        {
            scalar dalpha = dPI/(C1 + C2);

            const scalar p10 = p1;
            const scalar p20 = p2;
            p1 = p10 - dalpha*C1;
            p2 = p20 + dalpha*C2;
            dPI = p1 - p2;

            while
            (
                alpha1 + dalpha > 1
             || alpha1 + dalpha < 0
             || mag(dPI) > mag(dPIOld)
             || (p1 < small && p2 < small)
            )
            {
                dalpha /= 2.0;

                p1 = p10 - dalpha*C1;
                p2 = p20 + dalpha*C2;
                dPI = p1 - p2;
            }

            alpha1 += dalpha;
            alpha2 = 1.0 - alpha1;

            rho1 = alphaRho1/max(alpha1, rAlpha1);
            rho2 = alphaRho2/max(alpha2, rAlpha2);

            e1 = thermo1.calcCelle(p1, celli);
            e2 = thermo2.calcCelle(p2, celli);

            PI = fluid_.cellPI(celli);

            if (mag(dPI/PI) < pRelTol)
            {
                break;
            }

            C1 = sqr(phase1.cellSpeedOfSound(PI, celli))*rho1/alpha1;
            C2 = sqr(phase2.cellSpeedOfSound(PI, celli))*rho2/alpha2;

            dPIOld = dPI;
        } while (++iter < maxIter);

        phase1.alphaRhoE()[celli] = alpha1*rho1*(e1 + KE1);
        phase2.alphaRhoE()[celli] = alpha2*rho2*(e2 + KE2);
    }
}


void Foam::instantPressureRelaxation::FX
(
    const ScalarMultivariateEquation::VarType& alphaPI,
    const label li,
    scalarList& fx
) const
{
    fx.setSize(nEqns());
    fx = 0.0;

    const scalar PI = alphaPI.last();
    scalar sumAlpha = 0;

    forAll(rho0_, i)
    {
        const label phasei = phases_[i];

        const scalar alpha = alphaPI[i];
        const scalar rho = alphaRho0_[i]/max(alpha, phaseModels_[phasei].residualAlpha().value());

        // const scalar rho = alphaPI[i];
        // const scalar alpha = alphaRho0_[i]/max(rho, small);

        thermos_[phasei].rhoRef()[li] = rho;

        const scalar e = thermos_[phasei].calcCelle(PI, li);
        thermos_[phasei].he()[li] = e;

        fx[i] =
        (
            alphaRho0_[i]*(e - e0_[i])
          + 0.5*(PI + PI0_)*(alpha - alpha0_[i])
        );///max(PI, small);

        // fx[i] =
        //     (
        //         2.0*rho*rho0_[i]*(e - e0_[i])
        //       + (PI + PI0_)*(rho - rho0_[i])
        //     );

        sumAlpha += alpha;
    }
    fx.last() = (sumAlpha - alphaMax_);
}


// void Foam::instantPressureRelaxation::jacobian
// (
//     const typename ScalarMultivariateEquation::VarType& alphaPI,
//     const label li,
//     scalarList& fx,
//     RectangularMatrix<scalar>& J
// ) const
// {
//     fx.setSize(nEqns());
//     fx = 0.0;
//     J = Zero;
//
//     const scalar PI = alphaPI.last();
//     scalar sumAlpha = 0;
//
//     forAll(rho0_, i)
//     {
//         const label phasei = phases_[i];
//         const scalar alpha = alphaPI[i];
//         const scalar rho = alphaRho0_[i]/max(alpha, small);
//         thermos_[phasei].rhoRef()[li] = rho;
//
//         const scalar e = thermos_[phasei].calcCelle(PI, li);
//         thermos_[phasei].he()[li] = e;
//
//         fx[i] =
//             alphaRho0_[i]*(e - e0_[i])
//           + 0.5*(PI + PI0_)*(alpha - alpha0_[i]);
//         J(i, i) = 0.5*(PI + PI0_);
//         J(i, alphaPI.size()-1) = PI;
//
//         sumAlpha += alpha;
//     }
//     fx.last() = (sumAlpha - alphaMax_)*PI;
//     J(alphaPI.size()-1, alphaPI.size()-1) = sumAlpha - alphaMax_;
// }


bool Foam::instantPressureRelaxation::solve
(
    const scalar& deltaT
)
{
    if (!useRootSolver_)
    {
        relax4(deltaT);
        return true;
    }


    DynamicList<scalar> KEs(phaseModels_.size());
    DynamicList<scalar> alphaPI(phaseModels_.size() + 1);
    DynamicList<scalar> lower(phaseModels_.size() + 1, 0.0);
    DynamicList<scalar> upper(phaseModels_.size() + 1, great);
    DynamicList<scalar> dX(phaseModels_.size() + 1);
    const scalarField relTols(rootSolver_->relTolerances());
    const scalarField absTols(rootSolver_->absTolerances());
    DynamicList<label> inds(phaseModels_.size() + 1);
    forAll(phaseModels_[0], celli)
    {
        // Clear temporary variables
        alpha0_.clear();
        alphaRho0_.clear();
        rho0_.clear();
        e0_.clear();
        KEs.clear();
        alphaPI.clear();
        phases_.clear();
        dX.clear();
        lower.clear();
        upper.clear();
        inds.clear();

        // Look for non-negligible phases
        alphaMax_ = 1.0;
        forAll(fixedPhaseModels_, i)
        {
            alphaMax_ -= fixedPhaseModels_[i][celli];
        }

        forAll(phaseModels_, phasei)
        {
            const phaseModel& phase = phaseModels_[phasei];
            const scalar alpha = phase[celli];
            const scalar alphaRho = phase.alphaRho()[celli];
            if
            (
                alpha > phase.residualAlpha().value()
             && alphaRho > phase.residualAlphaRho().value()
            )
            {
                fluidBlastThermo& thermo = thermos_[phasei];
                const scalar rho = alphaRho/alpha;
                thermo.rhoRef()[celli] = rho;

                scalar& e = thermo.he()[celli];
                KEs.append(0.5*magSqr(phase.alphaRhoU()[celli]/alphaRho));

                e = phase.alphaRhoE()[celli]/alphaRho - KEs.last();

                scalar& T = thermo.T()[celli];
                scalar eTest = thermo.cellhe(thermo.TLow(), celli);
                if (e < eTest)
                {
                    e = eTest;
                    T = thermo.TLow();
                }
                else
                {
                    T = thermo.cellThe(e, T, celli);
                }

                thermo.p()[celli] = phase.cellpRhoT(celli);

                phases_.append(phasei);

                alpha0_.append(alpha);
                alphaRho0_.append(alphaRho);
                rho0_.append(rho);
                e0_.append(e);

                alphaPI.append(alpha);
                dX.append(phase.residualAlpha().value());
                lower.append(0.0);
                upper.append(1.0);
                inds.append(phasei);

                // alphaPI.append(rho);
                // dX.append(rho*1e-4);
                // lower.append(0.0);
                // upper.append(great);
                // inds.append(phasei);
            }
        }

        if (phases_.size() < 2) continue;

        // Update interfacial pressure
        PI0_ = fluid_.cellPI(celli);

        // Check if pressure actually needs to be relaxed
        scalar maxdPI = 0.0;
        forAll(phases_, i)
        {
            maxdPI = max
            (
                maxdPI,
                mag(thermos_[phases_[i]].p()[celli] - PI0_)
            );
        }

        if (maxdPI < PI0_*1e-6)
        {
            continue;
        }

        // Add interfacial pressure info
        alphaPI.append(PI0_);
        dX.append(max(PI0_*1e-4, 1.0));
        lower.append(small);
        upper.append(great);
        inds.append(phaseModels_.size());


        // Update number of variables, equations, and set bounds and dX terms
        pressureRelaxationSolver::nEqns_ = alphaPI.size();
        ScalarMultivariateEquation::nEqns_ = alphaPI.size();
        ScalarMultivariateEquation::nVar_ = ScalarMultivariateEquation::nEqns_;

        this->dX_ = dX;
        this->setLowerLimits(lower);
        this->setUpperLimits(upper);

        // Update tolerances (subset of original tolerances)
        rootSolver_->absTolerances() = scalarField(absTols, inds);
        rootSolver_->relTolerances() = scalarField(relTols, inds);

        // Root solve
        alphaPI = rootSolver_->solve(alphaPI, celli);

        // Transfer final information
        forAll(rho0_, i)
        {
            const label phasei = phases_[i];
            phaseModel& phase = phaseModels_[phasei];
            fluidBlastThermo& thermo = thermos_[phasei];

            const scalar alphaRho = alphaRho0_[i];

            phase[celli] = alphaPI[i];
            const scalar rho =
                alphaRho0_[i]/max(alphaPI[i], phase.residualAlpha().value());

            // const scalar rho = alphaPI[i];
            // phase[celli] = alphaRho0_[i]/max(rho, small);

            thermo.rhoRef()[celli] = rho;
            phase.he()[celli] =
                max
                (
                    thermo.calcCelle(alphaPI.last(), celli),
                    thermo.cellhe(thermo.TLow(), celli)
                );
            phase.alphaRhoE()[celli] = alphaRho*(phase.he()[celli] + KEs[i]);
        }
    }

    rootSolver_->absTolerances() = absTols;
    rootSolver_->relTolerances() = relTols;

    return true;
}

// void Foam::instantPressureRelaxation::FX
// (
//     const ScalarMultivariateEquation::VarType& rhoPI,
//     const label li,
//     scalarList& fx
// ) const
// {
//     fx.setSize(nEqns());
//     fx = 0.0;
//
//     const scalar PI = rhoPI.last();
//     scalar sumAlpha = 0;
//
//     forAll(rho0_, i)
//     {
//         const label phasei = phases_[i];
//         const scalar alpha = rhoPI[i];
//         const scalar rho = rhoPI[i + phases_.size()];
//         thermos_[phasei].rhoRef()[li] = rho;
//         const scalar e = thermos_[phasei].calcCelle(PI, li);
//
//         fx[i] = alpha*rho - alphaRho0_[i];
//         fx[i + phases_.size()] =
//             (e*alpha*rho - e0_[i]*alphaRho0_[i])
//           + (alpha - alpha0_[i])*0.5*(PI0_ + PI);
//
//         sumAlpha += alpha;
//     }
//     fx.last() = sumAlpha - alphaMax_;
// }
//
//
// bool Foam::instantPressureRelaxation::solve
// (
//     const scalar& deltaT
// )
// {
//     // // relax1(deltaT);
//     if (!useRootSolver_)
//     {
//         relax4(deltaT);
//         return true;
//     }
//
//
//     DynamicList<scalar> KEs(phaseModels_.size());
//     DynamicList<scalar> rhoPI(2*phaseModels_.size() + 1);
//     DynamicList<scalar> lower(2*phaseModels_.size() + 1, 0.0);
//     DynamicList<scalar> upper(2*phaseModels_.size() + 1, great);
//     DynamicList<scalar> dX(2*phaseModels_.size() + 1);
//     const scalarField relTols(rootSolver_->relTolerances());
//     const scalarField absTols(rootSolver_->absTolerances());
//     DynamicList<label> inds(2*phaseModels_.size());
//     forAll(phaseModels_[0], celli)
//     {
//         alpha0_.clear();
//         alphaRho0_.clear();
//         rho0_.clear();
//         e0_.clear();
//         KEs.clear();
//         rhoPI.clear();
//         phases_.clear();
//         dX.clear();
//         lower.clear();
//         upper.clear();
//         inds.clear();
//
//         // Look for non-negligible phases
//         forAll(phaseModels_, phasei)
//         {
//             const phaseModel& phase = phaseModels_[phasei];
//             const scalar alpha = phase[celli];
//             const scalar alphaRho = phase.alphaRho()[celli];
//             if
//             (
//                 alpha > phase.residualAlpha().value()
//              && alphaRho > phase.residualAlphaRho().value()
//             )
//             {
//                 fluidBlastThermo& thermo = thermos_[phasei];
//                 const scalar rho =
//                     alphaRho/max(alpha, phase.residualAlpha().value());
//                 thermo.rhoRef()[celli] = rho;
//
//                 scalar& e = thermo.he()[celli];
//                 KEs.append
//                 (
//                     0.5*magSqr(phase.alphaRhoU()[celli]/max(alphaRho, 1e-6))
//                 );
//
//                 e = phase.alphaRhoE()[celli]/max(alphaRho, 1e-6) - KEs.last();
//
//                 scalar& T = thermo.T()[celli];
//
//                 scalar eTest = thermo.cellhe(small, celli);
//                 if (e < eTest)
//                 {
//                     e = eTest;
//                     T = small;
//                 }
//                 else
//                 {
//                     T = thermo.cellThe(e, T, celli);
//                 }
//
//                 thermo.p()[celli] = phase.cellpRhoT(celli);
//
//                 alpha0_.append(alpha);
//                 alphaRho0_.append(alphaRho);
//                 rho0_.append(rho);
//                 e0_.append(e);
//                 rhoPI.append(alpha);
//
//                 phases_.append(phasei);
//                 inds.append(phasei);
//
//                 dX.append(max(alpha*1e-4, 1e-6));
//                 upper.append(1.0);
//                 // lower.append(phase.residualRho().value());
//             }
//         }
//         if (phases_.size() < 2) continue;
//
//         // Update interfacial pressure
//         PI0_ = fluid_.PI()[celli];
//
//         // Check if pressure actually needs to be relaxed
//         scalar maxdPI = 0.0;
//         forAll(phases_, i)
//         {
//             maxdPI = max
//             (
//                 maxdPI,
//                 mag(thermos_[phases_[i]].p()[celli] - PI0_)
//             );
//         }
//
//         if (maxdPI < PI0_*1e-6)
//         {
//             continue;
//         }
//
//         // Add the additional density information
//         forAll(phases_, i)
//         {
//             rhoPI.append(rho0_[i]);
//             dX.append(max(rho0_[i]*1e-4, phaseModels_[phases_[i]].residualRho().value()));
//             upper.append(great);
//             inds.append(phases_[i] + phaseModels_.size());
//         }
//
//         // Add interfacial pressure info
//         rhoPI.append(PI0_);
//         dX.append(PI0_*1e-4);
//         upper.append(great);
//         inds.append(2*phaseModels_.size());
//         // lower.append(0);
//
//         // Update number of variables, equations, and set bounds and dX terms
//         ScalarMultivariateEquation::nEqns_ = 2*phases_.size() + 1;
//         ScalarMultivariateEquation::nVar_ = ScalarMultivariateEquation::nEqns_;
//
//         lower.setSize(this->nVar());
//         // upper.setSize(this->nVar());
//
//         this->dX_ = dX;
//         this->setLowerLimits(lower);
//         this->setUpperLimits(upper);
//
//         // Update tolerances
//         rootSolver_->absTolerances() = scalarField(absTols, inds);
//         rootSolver_->relTolerances() = scalarField(relTols, inds);
//
//
//         rhoPI = rootSolver_->solve(rhoPI, celli);
//         forAll(rho0_, i)
//         {
//             const label phasei = phases_[i];
//             phaseModel& phase = phaseModels_[phasei];
//             fluidBlastThermo& thermo = thermos_[phasei];
//
//             const scalar alphaRho = alphaRho0_[i];
//             phase[celli] = rhoPI[i];
//             thermo.rhoRef()[celli] = rhoPI[i + phases_.size()];
//             phase.he()[celli] = thermo.calcCelle(rhoPI.last(), celli);
//             phase.alphaRhoE()[celli] =
//                 alphaRho*(phase.he()[celli] + KEs[i]);
//         }
//     }
//
//     rootSolver_->absTolerances() = absTols;
//     rootSolver_->relTolerances() = relTols;
//
//     // relax1(deltaT);
//     // relax4(deltaT);
//     return true;
// }
// ************************************************************************* //
