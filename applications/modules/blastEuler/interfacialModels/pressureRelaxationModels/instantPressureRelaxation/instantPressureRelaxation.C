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
        2*phaseModels_.size(),
        scalarField(2*phaseModels_.size(), 0.0),
        scalarField(2*phaseModels_.size(), great)
    ),
    alphaPIEqn_(phaseModels_, thermos_),
    simpleEqn_(phaseModels_, thermos_, alpha0_, e0_, Pi0_)
    // rootSolver_(rootSolver::New(*this, dict)),
    // uniRootSolver_(univariateRootSolver::New(simpleEqn_, dict)),
    // integrator_(scalarIntegrator::New(alphaPIEqn_, dict))
{
    scalarField ll(this->lowerLimits());
    for (label i = phaseModels_.size(); i < ll.size(); i++)
    {
        ll[i] = -great;
    }
    this->setLowerLimits(ll);

    scalarField ul(this->upperLimits());
    for (label i = 0; i < phaseModels_.size(); i++)
    {
        ul[i] = 1.0;
    }
    this->setUpperLimits(ul);

    pressureRelaxationSolver::nEqns_ = phaseModels_.size()*2;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::instantPressureRelaxation::~instantPressureRelaxation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::instantPressureRelaxation::FX
(
    const ScalarMultivariateEquation::VarType& alphaPI,
    const label li,
    scalarList& fx
) const
{
    // fx.setSize(nEqns());
    // fx = 0.0;
    // scalar sumAlpha = 0;
    //
    // forAll(phaseModels_, phasei)
    // {
    //     scalar alphaRho = phaseModels_[phasei].alphaRho()[li];
    //     if (alphaRho < 1e-10)
    //     {
    //         continue;
    //     }
    //     const scalar alpha = alphaPI[phasei];
    //     const scalar rho = alphaRho/max(alpha, 1e-6);
    //     thermos_[phasei].rhoRef()[li] = rho;
    //     scalar e = thermos_[phasei].calcCelle(alphaPI.last(), li);
    //     thermos_[phasei].he()[li] = e;
    //     const scalar pi = thermos_[phasei].cellpRhoT(li, false);
    //
    //     fx[phasei] = pi - alphaPI.last();
    //
    //         // 2.0*rhoPI[phasei]*rho0_[phasei]*(e - e0_[phasei])
    //       // + (rhoPI.last() + PI0_)*(rhoPI[phasei] - rho0_[phasei]);
    //     sumAlpha += alpha;
    // }
    // fx.last() = 1.0 - sumAlpha;
}


void Foam::instantPressureRelaxation::jacobian
(
    const ScalarMultivariateEquation::VarType& rhoPI,
    const label li,
    scalarList& fx,
    RectangularMatrix<scalar>& J
) const
{
    // J.setSize(nEqns(), nEqns());
    // fx.setSize(nEqns(), 0.0);
    // scalar sumAlpha = 0;
    // fx = 0.0;
    // J = Zero;
    //
    // forAll(phaseModels_, phasei)
    // {
    //     const scalar alphaRho = phaseModels_[phasei].alphaRho()[li];
    //     if (alphaRho < 1e-10)
    //     {
    //         continue;
    //     }
    //     const scalar rhos = max(rhoPI[phasei], 1e-10);
    //     const scalar alpha = alphaRho/rhos;
    //     thermos_[phasei].rhoRef()[li] = rhoPI[phasei];
    //     scalar e = thermos_[phasei].calcCelle(rhoPI.last(), li);
    //     thermos_[phasei].he()[li] = e;
    //
    //     const scalar pi = thermos_[phasei].cellpRhoT(li, false);
    //     const scalar dpdRho = thermos_[phasei].celldpdRho(li);
    //     const scalar dAlphadRho = -alphaRho/sqr(rhos);
    //
    //     fx[phasei] = pi - rhoPI.last();
    //     J(phasei, phasei) = thermos_[phasei].celldpdRho(li);
    //     J(phasei, phaseModels_.size()) = -1.0;
    //     J(phaseModels_.size(), phasei) = dAlphadRho;
    //
    //     sumAlpha += alpha;
    // }
    //
    // // forAll(phaseModels_, phasei)
    // // {
    // //     scalar alphaRho = phaseModels_[phasei].alphaRho()[li];
    // //     if (alphaRho < 1e-10)
    // //     {
    // //         continue;
    // //     }
    // //     scalar alpha = alphaRho/max(rhoPI[phasei], 1e-10);
    // //     thermos_[phasei].rho()[li] = rhoPI[phasei];
    // //     scalar e = thermos_[phasei].calcCelle(rhoPI.last(), li);
    // //     thermos_[phasei].he()[li] = e;
    // //
    // //     fx[phasei] =
    // //         2.0*rhoPI[phasei]*rho0_[phasei]*(e - e0_[phasei])
    // //       + (rhoPI.last() + PI0_)*(rhoPI[phasei] - rho0_[phasei]);
    // //
    // //     J(phasei, phasei) =
    // //         2.0*rho0_[phasei]
    // //        *(
    // //             (e - e0_[phasei])
    // //           + rhoPI[phasei]*thermos_[phasei].celldpde(li)
    // //         );
    // //     J(phasei, phaseModels_.size()) = -alphaRho/sqr(rhoPI[phasei]);
    // //     J(phaseModels_.size(), phasei) = rho0_[phasei] - rhoPI[phasei];
    // //
    // //     sumAlpha += alpha;
    // // }
    // // fx.last() = sumAlpha - 1.0;
    // J(phaseModels_.size(), phaseModels_.size()) = small;
}


void Foam::instantPressureRelaxation::relax1
(
    const scalar deltaT
)
{
    const bool print = false;
    scalar D = 2.0;
    scalar eps = 0.5;
    const label maxIter = 1000;
    const scalar pRelTol = 1e-8;
    const scalar pAbsTol = 1.0;
    phaseModel& phase1 = phaseModels_[0];
    phaseModel& phase2 = phaseModels_[1];

    fluidBlastThermo& thermo1 = thermos_[0];
    fluidBlastThermo& thermo2 = thermos_[1];

    const scalar rAlpha1 = phase1.residualAlpha().value();
    const scalar rAlpha2 = phase2.residualAlpha().value();

    bool end = false;

    //Procedure 1
    forAll(phase1, celli)
    {
        scalar& alpha1 = phase1[celli];
        scalar& alpha2 = phase2[celli];

        if (alpha1 < rAlpha1 || alpha2 < rAlpha2)
        {
            continue;
        }

        const scalar alpha10 = alpha1;
        const scalar alpha20 = alpha2;

        const scalar alphaRho1 = phase1.alphaRho()[celli];
        const scalar alphaRho2 = phase2.alphaRho()[celli];

        scalar& rho1 = thermo1.rhoRef()[celli];
        scalar& rho2 = thermo2.rhoRef()[celli];

        rho1 = alphaRho1/max(alpha1, rAlpha1);
        rho2 = alphaRho2/max(alpha2, rAlpha2);

        scalar& e1 = thermo1.he()[celli];
        scalar& e2 = thermo2.he()[celli];

        scalar& T1 = thermo1.T()[celli];
        scalar& T2 = thermo2.T()[celli];

        scalar& alphaRhoE1 = phase1.alphaRhoE()[celli];
        scalar& alphaRhoE2 = phase2.alphaRhoE()[celli];

        const scalar alphaRhoE10 = alphaRhoE1;
        const scalar alphaRhoE20 = alphaRhoE2;

        const scalar KE1 = 0.5*magSqr(phase1.alphaRhoU()[celli]/alphaRho1);
        const scalar KE2 = 0.5*magSqr(phase2.alphaRhoU()[celli]/alphaRho2);

        e1 = (alphaRhoE1)/max(alphaRho1, rAlpha1) - KE1;
        e2 = (alphaRhoE2)/max(alphaRho2, rAlpha2) - KE2;

        // T1 = max(thermo1.cellThe(e1, T1, rho1), small);
        // T2 = max(thermo2.cellThe(e2, T2, rho2), small);

        scalar& p1 = thermo1.p()[celli];
        scalar& p2 = thermo2.p()[celli];

        p1 = phase1.cellpRhoT(celli);
        p2 = phase2.cellpRhoT(celli);

        scalar PI = fluid_.cellPI(celli);

        scalar dPI = p1 - p2;

        if (mag(dPI)/PI < pRelTol)
        {
            continue;
        }

        scalar dPIOld = dPI;

        scalar dalpha =
            eps
           *min
            (
                max(rAlpha1, alpha1),
                max(rAlpha2, alpha2)
            );
        if (alpha1 + dalpha > 1)
        {
            dalpha *= -1.0;
        }
        const scalar dalpha0 = dalpha;

        scalar delta = 0.0;

        if(print)
            Info<<"alpha0: "<<alpha10<<" "<<alpha20<<endl
                <<"rho0: "<<rho1<<" "<<rho2<<endl
                <<"p0: "<<p1<<" "<<p2<<" "<<PI<<endl;

        if (p1 < small || p2 < small)
        {
            end = true;
        }
        label iter = 0;
        OStringStream os;
        do
        {
            alpha1 += dalpha;
            alpha2 = 1.0 - alpha1;

            delta += dalpha*PI;

            rho1 = alphaRho1/max(alpha1, rAlpha1);
            rho2 = alphaRho2/max(alpha2, rAlpha2);

            alphaRhoE1 = alphaRhoE10 - delta;
            alphaRhoE2 = alphaRhoE20 + delta;

            e1 = alphaRhoE1/max(alphaRho1, rAlpha1) - KE1;
            e2 = alphaRhoE2/max(alphaRho2, rAlpha2) - KE2;

            if(print)
                Info<<"iter " << iter<<":"<<nl
                << "    alpha: " << alpha1 << " "<<alpha2<<nl
                << "    alphaRho: " << alphaRho1 << " "<<alphaRho2<<nl
                << "    alphaRhoE: " << alphaRhoE1 << " "<<alphaRhoE2<<nl
                << "    alphaRhoe: " << alphaRhoE1-KE1 << " "<<alphaRhoE2-KE2<<nl
                << "    rho: " << rho1 << " "<<rho2<<nl
                << "    e: " << e1 << " "<<e2<<endl
                << "    T: " << T1 << " "<<T2<<endl;

            // T1 = max(thermo1.cellThe(e1, T1, celli), small);
            // T2 = max(thermo2.cellThe(e2, T2, celli), small);

            p1 = thermo1.cellpRhoT(celli);
            p2 = thermo2.cellpRhoT(celli);
            PI = fluid_.cellPI(celli);

            dPI = p1 - p2;

            if(print)Info<<"    p: " <<PI <<" "<<p1<<" "<<p2<<nl
                <<"    error: "<<dPI<<" "<<dPI/PI<<endl;

            if (mag(dPI/PI) < pRelTol || D < small)
            {
                break;
            }

            // Info<<D<<" "<<dPI<<" "<<dPIOld<<endl;
            if (dPI*dPIOld < 0)
            {
                dalpha = dalpha/D;
            }
            else if (mag(dPI) > mag(dPIOld))
            {
                dalpha = -dalpha/D;
            }
            dPIOld = dPI;

            // Info<<"   dalpha: "<<dalpha<<" "<<alpha1+dalpha<<endl;
            if (alpha1 + dalpha >= 1 - mag(dalpha))
            {
                dalpha = (1.0 - alpha1)/D;
                // Info<<"   dalpha+: "<<dalpha<<endl;
            }
            if (alpha1 + dalpha <= mag(dalpha))
            {
                dalpha = alpha1/D;
                // Info<<"   dalpha-: "<<dalpha<<endl;
            }
        } while (++iter < maxIter);

        phase1.alphaRhoE()[celli] = alpha1*rho1*(e1 + KE1);
        phase2.alphaRhoE()[celli] = alpha2*rho2*(e2 + KE2);
    }
}


void Foam::instantPressureRelaxation::relax4
(
    const scalar deltaT
)
{
    const bool print = false;
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

    bool end = false;

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
         // || alphaRho1 < rAlphaRho1
         // || alphaRho2 < rAlphaRho2
        )
        {
            continue;
        }

        OStringStream os;

        const scalar alpha10 = alpha1;
        const scalar alpha20 = alpha2;

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

        if(print)
            os<<"alpha0: "<<alpha10<<" "<<alpha20<<endl
                <<"rho0: "<<rho1<<" "<<rho2<<endl
                <<"alphaRho0: "<<alphaRho1<<" "<<alphaRho2<<endl
                <<"e0: "<<e1<<" "<<e2<<endl;

        scalar e1Test = thermo1.cellhe(small, celli);
        if (e1 < e1Test)
        {
            e1 = e1Test;
            T1 = small;
            if (print) Info<<"corr1: "<<e1<<endl;
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
            if (print) Info<<"corr2: "<<e2<<endl;
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

        if(print)
            Pout<<"p0: "<<PI<<" "<<p1<<" "<<p2<<endl
                <<"e0: "<<e1<<" "<<e2<<endl
                <<"T0: "<<T1<<" "<<T2<<endl;

        label iter = 0;
        label iter1 = 0;
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

            if(print)
                Pout<<"iter " << iter1<<":"<<nl
                << "    alpha: " << alpha1 << " "<<alpha2<<nl
                << "    alphaRho: " << alphaRho1 << " "<<alphaRho2<<nl
                << "    rho: " << rho1 << " "<<rho2<<nl
                << "    p: " <<PI <<" "<<p1<<" "<<p2<<nl
                << "    e: " << e1 << " "<<e2<<endl
                << "    T: " << T1 << " "<<T2<<endl;



            e1 = thermo1.calcCelle(p1, celli);
            e2 = thermo2.calcCelle(p2, celli);

            // T1 = max(thermo1.cellThe(e1, T1, celli), small);
            // T2 = max(thermo2.cellThe(e2, T2, celli), small);

            // p1 = thermo1.cellpRhoT(celli);
            // p2 = thermo2.cellpRhoT(celli);

            PI = fluid_.cellPI(celli);

            const scalar dPIRel = dPI/PI;

            if(print)
                os<<"    p: " <<PI <<" "<<p1<<" "<<p2<<nl
                <<"    error: "<<dPI<<" "<<dPIRel<<endl;

            if (mag(dPIRel) < pRelTol)
            {
                break;
            }

            C1 = sqr(phase1.cellSpeedOfSound(PI, celli))*rho1/alpha1;
            C2 = sqr(phase2.cellSpeedOfSound(PI, celli))*rho2/alpha2;
            iter1++;

            dPIOld = dPI;
        } while (++iter < maxIter);

        phase1.alphaRhoE()[celli] = alpha1*rho1*(e1 + KE1);
        phase2.alphaRhoE()[celli] = alpha2*rho2*(e2 + KE2);

        if (mag(dPI) > small && print && iter1 > 1)
        {
            os<<"nIter = " << iter<<": "<<alpha1<<" "<<alpha10<<endl;
            Info<<word(os.str())<<endl;
            end = true;
        }
    }
    if (end)
        std::exit(0);
}


bool Foam::instantPressureRelaxation::solve
(
    const scalar& deltaT
)
{
    // alpha0_.setSize(phaseModels_.size());
    // e0_.setSize(phaseModels_.size());
    // Pi0_.setSize(phaseModels_.size());
    //
    // // scalarField alpha_e(phaseModels_.size() + 1);
    // forAll(phaseModels_[0], celli)
    // {
    //     // if (celli !=84)
    //     // {
    //     //     continue;
    //     // }
    //     forAll(phaseModels_, phasei)
    //     {
    //         alpha0_[phasei] = phaseModels_[phasei][celli];
    //         e0_[phasei] = phaseModels_[phasei].he()[celli];
    //         Pi0_[phasei] = thermos_[phasei].p()[celli];
    //         // alphaPI[phasei] = alpha0_[phasei];
    //     }
    //     // PI0_ = fluid_.PI()[celli];
    //     // alphaPI.last() = fluid;
    //
    //     // alphaPI = rootSolver_->solve(alphaPI, celli);
    //
    //     scalar alpha = alpha0_[0];
    //     // if (alpha0_[0] > 1e-6 && alpha0_[0] < 1.0 - 1e-6)
    //     {
    //         alpha = uniRootSolver_->solve(alpha, 0.0, 1.0, celli);
    //     }
    // }
    // // std::exit(0);

    // relax1(deltaT);
    relax4(deltaT);
    return true;
}

// ************************************************************************* //
