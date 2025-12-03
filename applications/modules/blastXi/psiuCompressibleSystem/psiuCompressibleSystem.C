/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2025
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "psiuCompressibleSystem.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(psiuCompressibleSystem, 0);
    addToRunTimeSelectionTable
    (
        compressibleSystem,
        psiuCompressibleSystem,
        singlePhase
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::psiuCompressibleSystem::psiuCompressibleSystem
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    compressibleSystem(dict, mesh),
    thermo_(psiuMulticomponentThermo::New(mesh)),
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
        thermo_->rho()
    ),
    p_(thermo_->p()),
    T_(thermo_->T()),
    e_(thermo_->he()),
    eu_(thermo_->heu()),
    rhoEu_
    (
        IOobject
        (
            "rhoEu",
            mesh.time().name(),
            mesh
        ),
        eu_*rho_
    ),

    combustionProperties_
    (
        IOobject
        (
            "combustionProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    b_(thermo_->Y("b")),
    Xi_
    (
        IOobject
        (
            "Xi",
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    Su_
    (
        IOobject
        (
            "Su",
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    St_
    (
        IOobject
        (
            "St",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Xi_*Su_
    ),

    SuMin_(0.01*Su_.average()),
    SuMax_(4.0*Su_.average()),

    unstrainedLaminarFlameSpeed_(laminarFlameSpeed::New(thermo_())),

    SuModel_(combustionProperties_.lookup("SuModel")),
    XiModel_(combustionProperties_.lookup("XiModel")),

    sigmaExt_(combustionProperties_.lookup("sigmaExt")),
    XiCoef_(combustionProperties_.lookup("XiCoef")),
    XiShapeCoef_(combustionProperties_.lookup("XiShapeCoef")),
    uPrimeCoef_(combustionProperties_.lookup("uPrimeCoef")),

    ign_(combustionProperties_, mesh.time(), mesh)
{
    thermo_->validate("psiuCompressibleSystem", "ea");

    turbulence_ =
        compressible::momentumTransportModel::New
        (
            rho_,
            U_,
            rhoPhi_,
            thermo_()
        );
    thermophysicalTransport_.set
    (
        new turbulenceThermophysicalTransportModels::unityLewisEddyDiffusivity
        <
            RASThermophysicalTransportModel
            <
                ThermophysicalTransportModel
                <
                    compressibleMomentumTransportModel,
                    fluidThermo
                >
            >
        >
        (
            turbulence_(),
            thermo_(),
            true
        )
    );

    fluxScheme_ = fluxScheme::NewSingle(phi_);
    encode();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::psiuCompressibleSystem::~psiuCompressibleSystem()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar Foam::psiuCompressibleSystem::calcStCorr
(
    const volScalarField& c,
    const surfaceScalarField& nf,
    const dimensionedScalar& dMgb
) const
{
    dimensionedScalar StCorr("StCorr", dimless, 1.0);

    if (ign_.igniting())
    {
        // Calculate volume of ignition kernel
        const dimensionedScalar Vk
        (
            "Vk",
            dimVolume,
            gSum(c*mesh().V().primitiveField())
        );
        dimensionedScalar Ak("Ak", dimArea, 0.0);

        if (Vk.value() > small)
        {
            // Calculate kernel area from its volume
            // and the dimensionality of the case

            switch(mesh().nGeometricD())
            {
                case 3:
                {
                    // Assume it is part-spherical
                    const scalar sphereFraction
                    (
                        combustionProperties_.lookup<scalar>
                        (
                            "ignitionSphereFraction"
                        )
                    );

                    Ak = sphereFraction*4.0*constant::mathematical::pi
                       *pow
                        (
                            3.0*Vk
                           /(sphereFraction*4.0*constant::mathematical::pi),
                            2.0/3.0
                        );
                }
                break;

                case 2:
                {
                    // Assume it is part-circular
                    const dimensionedScalar thickness
                    (
                        combustionProperties_.lookup("ignitionThickness")
                    );

                    const scalar circleFraction
                    (
                        combustionProperties_.lookup<scalar>
                        (
                            "ignitionCircleFraction"
                        )
                    );

                    Ak = circleFraction*constant::mathematical::pi*thickness
                       *sqrt
                        (
                            4.0*Vk
                           /(
                               circleFraction
                              *thickness
                              *constant::mathematical::pi
                            )
                        );
                }
                break;

                case 1:
                    // Assume it is plane or two planes
                    Ak = dimensionedScalar
                    (
                        combustionProperties_.lookup("ignitionKernelArea")
                    );
                break;
            }

            // Calculate kernel area from b field consistent with the
            // discretisation of the b equation.
            const volScalarField mgb
            (
                fvc::div(nf, b_, "div(phiSt,b)") - b_*fvc::div(nf) + dMgb
            );
            const dimensionedScalar AkEst = gSum(mgb*mesh().V().primitiveField());

            StCorr.value() = max(min((Ak/AkEst).value(), 10.0), 1.0);

            Info<< "StCorr = " << StCorr.value() << endl;
        }
    }

    return StCorr;
}


void Foam::psiuCompressibleSystem::solve()
{
    turbulence_->predict();
    thermophysicalTransport_->predict();

    volSymmTensorField devTau(turbulence_->devTau());
    volScalarField divSigmaDotU
    (
        fvc::div
        (
            fvc::dotInterpolate(rho_.mesh().Sf(), turbulence_->devTau())
          & fluxScheme_->Uf()
        )
    );

    volScalarField deltaRho(fvc::div(rhoPhi_));
    volVectorField deltaRhoU
    (
        "deltaRhoU",
        fvc::div(rhoUPhi_)
      - g_*rho_
      + fvc::div(devTau)
    );
    volScalarField deltaRhoE
    (
        fvc::div(rhoEPhi_)
      - (rhoU_ & g_)
      - (thermophysicalTransport_->divq(e_) & e_)
      - divSigmaDotU
    );
    volScalarField deltaRhoEu
    (
        fvc::div(fluxScheme_->energyFlux(rho_, U_, eu_, p_))
      - (rhoU_ & g_)
      + (thermophysicalTransport_->divq(eu_) & eu_)
      - divSigmaDotU
    );

    // if (explicitViscosity_ && turbulence_.valid())
    {
        volSymmTensorField devTau(turbulence_->devTau());
        volScalarField divTauDotU
        (
            fvc::div
            (
                fvc::dotInterpolate(rho_.mesh().Sf(), devTau)
              & fluxScheme_->Uf()
            )
        );

        deltaRhoU += fvc::div(devTau);
        deltaRhoE +=
           (thermophysicalTransport_->divq(e_) & e_) + divTauDotU;
        deltaRhoEu +=
            (thermophysicalTransport_->divq(eu_) & eu_) + divTauDotU;
    }


    //- Store changed in mass, momentum and energy
    this->storeAndBlendDelta(deltaRho);
    this->storeAndBlendDelta(deltaRhoU);
    this->storeAndBlendDelta(deltaRhoE);
    this->storeAndBlendDelta(deltaRhoEu);

    //- Store old values
    const volScalarField rho0(rho_);
    this->storeAndBlendOld(rho_);

    this->storeAndBlendOld(rhoU_);
    this->storeAndBlendOld(rhoE_);
    this->storeAndBlendOld(rhoEu_);

    dimensionedScalar dT = rho_.time().deltaT();
    rho_ -= dT*deltaRho;
    rho_.correctBoundaryConditions();

    vector solutionDs((vector(rho_.mesh().solutionD()) + vector::one)/2.0);
    rhoU_ -= cmptMultiply(dT*deltaRhoU, solutionDs);
    rhoE_ -= dT*deltaRhoE;
    rhoEu_ -= dT*deltaRhoEu;

    if (thermo_->containsSpecie("ft"))
    {
        volScalarField& ft = thermo_->Y("ft");
        volScalarField deltaRhoFt
        (
            fvc::div(fluxScheme_->flux(ft, rhoPhi_))
        );

        volScalarField rhoft(rho_.name() + ft.name(), rho0*ft);
        this->storeAndBlendOld(rhoft);
        this->storeAndBlendDelta(deltaRhoFt);

        ft = (rhoft - dT*deltaRhoFt)/rho_;
        ft.max(0.0);
        ft.correctBoundaryConditions();
    }

    if (ign_.ignited())
    {
        const fvMesh& mesh = this->mesh();

        // progress variable
        // ~~~~~~~~~~~~~~~~~
        volScalarField c("c", scalar(1.0) - b_);

        // Unburnt gas density
        // ~~~~~~~~~~~~~~~~~~~
        volScalarField rhou(thermo_->rhou());

        // Calculate flame normal etc.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        volVectorField n("n", fvc::grad(b_));

        volScalarField mgb(mag(n));

        dimensionedScalar dMgb =
            1.0e-3
           *(b_*c*mgb)().weightedAverage(mesh.V())
           /((b_*c)().weightedAverage(mesh.V()) + small)
          + dimensionedScalar(mgb.dimensions(), small);
        mgb += dMgb;

        surfaceVectorField SfHat(mesh.Sf()/mesh.magSf());
        surfaceVectorField nfVec(fvc::interpolate(n));
        nfVec += SfHat*(fvc::snGrad(b_) - (SfHat & nfVec));
        nfVec /= (mag(nfVec) + dMgb);
        surfaceScalarField nf((mesh.Sf() & nfVec));
        n /= mgb;

        const volScalarField nDotn(n & n);

        const dimensionedScalar StCorr(calcStCorr(c, nf, dMgb));

        // Calculate turbulent flame speed flux
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        surfaceScalarField phiSt
        (
            "phiSt",
            fvc::interpolate(rhou*StCorr*Su_*Xi_)*nf
        );

        volScalarField deltaRhoB
        (
            "deltaRhoB",
            fvc::div((rhoPhi_ + phiSt)*fluxScheme_->interpolate(b_, "b"))
          - b_*fvc::div(phiSt)
        );

        volScalarField rhob(rho_.name() + b_.name(), rho0*b_);
        this->storeAndBlendOld(rhob);
        this->storeAndBlendDelta(deltaRhoB);
        b_ = (rhob - dT*deltaRhoB)/rho_;
        b_.maxMin(0.0, 1.0);
        b_.correctBoundaryConditions();


        tmp<volScalarField> tsigmas;
        if
        (
            SuModel_ == "equilibrium"
         || SuModel_ == "transport"
         || XiModel_ == "transport"
        )
        (
            tsigmas =
                (nDotn*fvc::div(phi_) - (n & fvc::grad(U_) & n))/Xi_
              + (
                    nDotn*fvc::div(Su_*n)
                  - (n & fvc::grad(Su_*n) & n)
                )*(Xi_ + scalar(1))/(2*Xi_)
        );

        // Calculate Xi flux
        // ~~~~~~~~~~~~~~~~~
        surfaceScalarField phiXi
        (
            phiSt
          - fvc::interpolate
            (
                fvc::laplacian(thermophysicalTransport_->DEff(b_), b_)/mgb
            )*nf
          + fvc::interpolate(rho_)*fvc::interpolate(Su_*(1.0/Xi_ - Xi_))*nf
        );

        const volScalarField SuOld(Su_);

        volScalarField Su0(unstrainedLaminarFlameSpeed_()());
        if (SuModel_ == "unstrained")
        {
            Su_ == Su0;
        }
        else if (SuModel_ == "equilibrium")
        {
            Su_ == Su0*max(scalar(1) - tsigmas/sigmaExt_, scalar(0.01));
        }
        else if (SuModel_ == "transport")
        {
            // Solve for the strained laminar flame speed
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            volScalarField deltaRhoSu
            (
                "deltaRhoSu",
                fvc::div((phi_ + phiXi)*fluxScheme_->interpolate(Su_, "Su"))
              - Su_*fvc::div(phiXi)
            );

            volScalarField rhoSu(rho_.name() + Su_.name(), rho0*Su_);
            this->storeAndBlendOld(rhoSu);
            this->storeAndBlendDelta(deltaRhoSu);


            Su_ = (rhoSu - dT*deltaRhoSu)/rho_;
            Su_.maxMin(SuMin_, SuMin_);
            Su_.correctBoundaryConditions();
        }
        else
        {
            FatalError
                << "Unknown Su model " << SuModel_
                << abort(FatalError);
        }


        // Calculate Xi according to the selected flame wrinkling model
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (XiModel_ == "fixed")
        {
            // Do nothing, Xi is fixed!
        }
        else if (XiModel_ == "algebraic")
        {
            volScalarField epsilon(pow(uPrimeCoef_, 3)*turbulence().epsilon());
            volScalarField tauEta(sqrt(thermo_->muu()/(rhou*epsilon)));

            volScalarField up(uPrimeCoef_*sqrt((2.0/3.0)*turbulence().k()));

            tmp<volScalarField> Reta
            (
                up
               /(
                    sqrt(epsilon*tauEta)
                  + dimensionedScalar(up.dimensions(), 1e-8)
                )
            );

            // Simple algebraic model for Xi based on Gulders correlation
            // with a linear correction function to give a plausible profile
            // for Xi
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            Xi_ ==
                1.0
                + (1.0 + (2*XiShapeCoef_)*(0.5 - b_))
                *XiCoef_*sqrt(up/(SuOld + SuMin_))*Reta;
        }
        else if (XiModel_ == "transport")
        {
            // Solve for the flame wrinkling
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            volScalarField deltaRhoXi
            (
                "deltaRhoXi",
                fvc::div
                (
                    (phi_ + phiXi)
                    *fluxScheme_->interpolate(Xi_, "Xi")
                )
                - fvc::div(phiXi)*Xi_
            );
            this->storeAndBlendDelta(deltaRhoXi);

            volScalarField rhoXi(rho_.name() + Xi_.name(), rho0*Xi_);
            this->storeAndBlendOld(Xi_);
            Xi_ = (rhoXi - dT*deltaRhoXi)/rho_;
            Xi_.max(1.0);
            Xi_.correctBoundaryConditions();
        }
        else
        {
            FatalError
                << "Unknown Xi model " << XiModel_
                << abort(FatalError);
        }
        St_ = Xi_*Su_;
    }
}


void Foam::psiuCompressibleSystem::storeExplicit()
{
    rhoAdvection_ = fvc::ddt(rho_);
    rhoUAdvection_ = fvc::ddt(rhoU_);
    rhoEAdvection_ = fvc::ddt(rhoE_);
    rhoEuAdvection_ = fvc::ddt(rhoEu_);

    if (thermo_->containsSpecie("ft"))
    {
        rhoftAdvection_ = fvc::ddt(rho_, thermo_->Y("ft"));
    }
    rhobAdvection_ = fvc::ddt(rho_, b_);

    if (SuModel_ == "transport")
    {
        rhoSuAdvection_ = fvc::ddt(rho_, Su_);
    }
    if (XiModel_ == "transport")
    {
        rhoXiAdvection_ = fvc::ddt(rho_, Xi_);
    }
}


void Foam::psiuCompressibleSystem::solveImplicit()
{
    if (turbulence_.valid())
    {
        turbulence_->predict();
    }
    if (thermophysicalTransport_.valid())
    {
        thermophysicalTransport_->predict();
    }

    tmp<surfaceVectorField> devTau;
    {
        tmp<fvVectorMatrix> divDevTau;
        if (!explicitViscosity_ && turbulence_.valid())
        {
            divDevTau =
                turbulence_->divDevTau(U_);
                // + fvc::grad((2.0/3.0)*rhoEff()*turbulence_->k());
        }

        // Solve momentum
        fvVectorMatrix UEqn
        (
            fvm::ddt(rhoEff(), U_)
          - rhoUAdvection_() // Change from advection
        ==
            models().source(rhoEff(), U_)
        );

        if (divDevTau.valid())
        {
            UEqn += divDevTau();
        }
        addUSource(UEqn);

        UEqn.relax();

        constraints().constrain(UEqn);
        UEqn.solve();
        constraints().constrain(U_);

        if (divDevTau.valid())
        {
            devTau = divDevTau().flux();
        }

        // Update kinetic energy and momentum
        K_ = 0.5*magSqr(U_);
        rhoU_ = rhoEff()*U_;
    }

    // Solve thermal energy diffusion
    {
        volScalarField& he = thermo().he();
        fvScalarMatrix EEqn
        (
            fvm::ddt(rhoEff(), he)
          - rhoEAdvection_()        // Explicit advection contribtion
          + fvc::ddt(rhoEff(), K_)  // Change in kinetic energy
         ==
            models().source(rhoEff(), he)
        );

        if (devTau.valid())
        {
            EEqn +=
                fvc::div(devTau & flux().Uf())
              + thermophysicalTransport_->divq(he);
        }

        EEqn.relax();

        constraints().constrain(EEqn);
        EEqn.solve();
        constraints().constrain(he);

        // Update total energy
        rhoE_ = rhoEff()*(he + K_);
    }

    // Solve unburnt thermal energy diffusion
    {
        fvScalarMatrix EuEqn
        (
            fvm::ddt(rhoEff(), eu_)
          - rhoEuAdvection_()        // Explicit advection contribtion
          + fvc::ddt(rhoEff(), K_)  // Change in kinetic energy
         ==
            models().source(rhoEff(), eu_)
        );

        if (devTau.valid())
        {
            EuEqn +=
                fvc::div(devTau & flux().Uf())
              + thermophysicalTransport_->divq(eu_);
        }

        EuEqn.relax();

        constraints().constrain(EuEqn);
        EuEqn.solve();
        constraints().constrain(eu_);

        // Update total energy
        rhoEu_ = rhoEff()*(eu_ + K_);
    }

    if (turbulence_.valid())
    {
        turbulence_->correct();
    }

    if (thermophysicalTransport_.valid())
    {
        thermophysicalTransport_->correct();
    }

    if (thermo_->containsSpecie("ft"))
    {
        volScalarField& ft = thermo_->Y("ft");
        fvScalarMatrix ftEqn
        (
            fvm::ddt(rhoEff(), ft) - rhoftAdvection_()
          + thermophysicalTransport_->divq(ft)
         ==
            models().source(rhoEff(), ft)
        );

        constraints().constrain(ftEqn);
        ftEqn.solve();
        constraints().constrain(ft);
    }

    if (ign_.ignited())
    {
        const fvMesh& mesh = this->mesh();

        // progress variable
        // ~~~~~~~~~~~~~~~~~
        volScalarField c("c", scalar(1.0) - b_);

        // Unburnt gas density
        // ~~~~~~~~~~~~~~~~~~~
        volScalarField rhou(thermo_->rhou());

        // Calculate flame normal etc.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        volVectorField n("n", fvc::grad(b_));

        volScalarField mgb(mag(n));

        dimensionedScalar dMgb =
            1.0e-3
           *(b_*c*mgb)().weightedAverage(mesh.V())
           /((b_*c)().weightedAverage(mesh.V()) + small)
          + dimensionedScalar(mgb.dimensions(), small);
        mgb += dMgb;

        surfaceVectorField SfHat(mesh.Sf()/mesh.magSf());
        surfaceVectorField nfVec(fvc::interpolate(n));
        nfVec += SfHat*(fvc::snGrad(b_) - (SfHat & nfVec));
        nfVec /= (mag(nfVec) + dMgb);
        surfaceScalarField nf((mesh.Sf() & nfVec));
        n /= mgb;

        const volScalarField nDotn(n & n);

        const dimensionedScalar StCorr(calcStCorr(c, nf, dMgb));

        // Calculate turbulent flame speed flux
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        surfaceScalarField phiSt
        (
            "phiSt",
            fvc::interpolate(rhou*StCorr*Su_*Xi_)*nf
        );

        scalar StCoNum = max
        (
            mesh.surfaceInterpolation::deltaCoeffs()
           *mag(phiSt)/(fvc::interpolate(rho_)*mesh.magSf())
        ).value()*rho_.time().deltaTValue();

        Info<< "Max St-Courant Number = " << StCoNum << endl;
        fvScalarMatrix bEqn
        (
            fvm::ddt(rho_, b_) - rhobAdvection_()
          - fvm::laplacian(thermophysicalTransport_->DEff(b_), b_)
        );

        forAll(ign_.sites(), i)
        {
            const ignitionSite& ignSite = ign_.sites()[i];

            if (ignSite.igniting())
            {
                forAll(ignSite.cells(), icelli)
                {
                    label ignCell = ignSite.cells()[icelli];
                    DebugInfo<< "Igniting cell " << ignCell;

                    DebugInfo<< " state :"
                        << ' ' << b_[ignCell]
                        << ' ' << Xi_[ignCell]
                        << ' ' << Su_[ignCell]
                        << ' ' << mgb[ignCell]
                        << endl;

                    bEqn.diag()[ignSite.cells()[icelli]] +=
                    (
                        ignSite.strength()*ignSite.cellVolumes()[icelli]
                       *rhou[ignSite.cells()[icelli]]/ignSite.duration()
                    )/(b_[ignSite.cells()[icelli]] + 0.001);
                }
            }
        }

        // Solve for b
        // ~~~~~~~~~~~
        bEqn.relax();

        constraints().constrain(bEqn);

        bEqn.solve();

        constraints().constrain(b_);

        Info<< "min(b) = " << min(b_).value() << endl;


        // Calculate Xi flux
        // ~~~~~~~~~~~~~~~~~
        tmp<surfaceScalarField> tphiXi;
        if (SuModel_ == "transport" || XiModel_ == "transport")
        {
            tphiXi =
                phiSt
              - fvc::interpolate
                (
                    fvc::laplacian(thermophysicalTransport_->DEff(b_), b_)/mgb
                )*nf
              + fvc::interpolate(rho_)*fvc::interpolate(Su_*(1.0/Xi_ - Xi_))*nf;
        }

        tmp<volScalarField> tsigmas;
        if
        (
            SuModel_ == "equilibrium"
         || SuModel_ == "transport"
         || XiModel_ == "transport"
        )
        (
            tsigmas =
                (nDotn*fvc::div(phi_) - (n & fvc::grad(U_) & n))/Xi_
              + (
                    nDotn*fvc::div(Su_*n)
                  - (n & fvc::grad(Su_*n) & n)
                )*(Xi_ + scalar(1))/(2*Xi_)
        );

        // Calculate Xi flux
        // ~~~~~~~~~~~~~~~~~
        surfaceScalarField phiXi
        (
            phiSt
          - fvc::interpolate
            (
                fvc::laplacian(thermophysicalTransport_->DEff(b_), b_)/mgb
            )*nf
          + fvc::interpolate(rho_)*fvc::interpolate(Su_*(1.0/Xi_ - Xi_))*nf
        );

        const volScalarField SuOld(Su_);

        volScalarField Su0(unstrainedLaminarFlameSpeed_()());
        if (SuModel_ == "unstrained")
        {
            Su_ == Su0;
        }
        else if (SuModel_ == "equilibrium")
        {
            Su_ == Su0*max(scalar(1) - tsigmas/sigmaExt_, scalar(0.01));
        }
        else if (SuModel_ == "transport")
        {
            const volScalarField& sigmas = tsigmas();
            volScalarField SuInf
            (
                Su0*max(scalar(1) - sigmas/sigmaExt_, scalar(0.01))
            );

            // Solve for the strained laminar flame speed
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            volScalarField Rc
            (
                (sigmas*SuInf*(Su0 - SuInf) + sqr(SuMin_)*sigmaExt_)
                /(sqr(Su0 - SuInf) + sqr(SuMin_))
            );

            fvScalarMatrix SuEqn
            (
                fvm::ddt(rho_, Su_) - rhoSuAdvection_()
            ==
              - fvm::SuSp(-rho_*Rc*Su0/Su_, Su_)
              - fvm::SuSp(rho_*(sigmas + Rc), Su_)
              + models().source(rho_, Su_)
            );

            SuEqn.relax();

            constraints().constrain(SuEqn);

            SuEqn.solve();

            constraints().constrain(Su_);

            // Limit the maximum Su
            // ~~~~~~~~~~~~~~~~~~~~
            Su_.min(SuMax_);
            Su_.max(SuMin_);
        }
        else
        {
            FatalError
                << "Unknown Su model " << SuModel_
                << abort(FatalError);
        }


        // Calculate Xi according to the selected flame wrinkling model
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (XiModel_ == "fixed")
        {
            // Do nothing, Xi is fixed!
        }
        else if (XiModel_ == "algebraic" || XiModel_ == "transport")
        {
            volScalarField epsilon(pow(uPrimeCoef_, 3)*turbulence().epsilon());
            volScalarField tauEta(sqrt(thermo_->muu()/(rhou*epsilon)));

            volScalarField up(uPrimeCoef_*sqrt((2.0/3.0)*turbulence().k()));
//             volScalarField up(sqrt(mag(diag(n * n) & diag(turbulence->r()))));

            tmp<volScalarField> Reta
            (
                up/(sqrt(epsilon*tauEta)+ dimensionedScalar(up.dimensions(), 1e-8))
            );
//             volScalarField l = 0.337*k*sqrt(k)/epsilon;
//             Reta *= max((l - dimensionedScalar(dimLength, 1.5e-3))/l, 0);

            if (XiModel_ == "algebraic")
            {
                // Simple algebraic model for Xi based on Gulders correlation
                // with a linear correction function to give a plausible profile
                // for Xi
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                Xi_ ==
                    1.0
                  + (1.0 + (2*XiShapeCoef_)*(0.5 - b_))
                   *XiCoef_*sqrt(up/(SuOld + SuMin_))*Reta;
            }
            else
            {
                volVectorField Ut(U_ + SuOld*Xi_*n);
                volScalarField sigmat
                (
                    nDotn*fvc::div(Ut) - (n & fvc::grad(Ut) & n)
                );

                // Calculate Xi transport coefficients based on Gulders correlation
                // and DNS data for the rate of generation
                // with a linear correction function to give a plausible profile for Xi
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                volScalarField XiEqStar
                (
                    1.001 + XiCoef_*sqrt(up/(SuOld + SuMin_))*Reta
                );
                volScalarField XiEq
                (
                    1.001
                 + (1.0 + (2*XiShapeCoef_)*(0.5 - b_))*(XiEqStar - 1.001)
                );

                volScalarField Gstar(0.28/tauEta);
                volScalarField R(Gstar*XiEqStar/(XiEqStar - 1.0));
                volScalarField G(R*(XiEq - 1.001)/XiEq);

                // R *= (Gstar + 2*mag(dev(symm(fvc::grad(U)))))/Gstar;

                // Solve for the flame wrinkling
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                fvScalarMatrix XiEqn
                (
                    fvm::ddt(rho_, Xi_) - rhoXiAdvection_()
                 ==
                    rho_*R
                  - fvm::Sp(rho_*(R - G), Xi_)
                  - fvm::Sp
                    (
                        rho_*max
                        (
                            sigmat - tsigmas,
                            dimensionedScalar(sigmat.dimensions(), 0)
                        ),
                        Xi_
                    )
                  + models().source(rho_, Xi_)
                );

                XiEqn.relax();

                constraints().constrain(XiEqn);

                XiEqn.solve();

                constraints().constrain(Xi_);

                // Correct boundedness of Xi
                // ~~~~~~~~~~~~~~~~~~~~~~~~~
                Xi_.max(1.0);
                Info<< "max(Xi) = " << max(Xi_).value() << endl;
                Info<< "max(XiEq) = " << max(XiEq).value() << endl;
            }
        }
        else
        {
            FatalError
                << "Unknown Xi model " << XiModel_
                << abort(FatalError);
        }
        St_ = Xi_*Su_;
    }

    this->decode();
    constraints().constrain(p_);
    p_.correctBoundaryConditions();

    turbulence_->correct();
    thermophysicalTransport_->correct();
}


void Foam::psiuCompressibleSystem::update()
{
    fluxScheme_->update
    (
        rho_,
        U_,
        e_,
        p_,
        speedOfSound()(),
        phi_,
        rhoPhi_,
        rhoUPhi_,
        rhoEPhi_
    );
}


void Foam::psiuCompressibleSystem::decode()
{
    U_.internalFieldRef() = rhoU_()/rho_();
    U_.correctBoundaryConditions();

    K_ = 0.5*magSqr(U_);

    e_.internalFieldRef() = rhoE_()/rho_() - K_();
    e_.correctBoundaryConditions();


    eu_.internalFieldRef() = rhoEu_()/rho_() - K_();
    forAll(b_, i)
    {
        if (b_[i] < 0.0001)
        {
            eu_[i] = e_[i];
        }
    }
    eu_.correctBoundaryConditions();

    thermo_->correct();
    p_.internalFieldRef() = rho_/thermo_->psi();
    p_.correctBoundaryConditions();
    rho_.boundaryFieldRef() ==
        thermo_->psi().boundaryField()*p_.boundaryField();

    rhoU_.boundaryFieldRef() = rho_.boundaryField()*U_.boundaryField();
    rhoE_.boundaryFieldRef() =
        rho_.boundaryField()*(e_.boundaryField() + K_.boundaryField());

    rhoEu_.boundaryFieldRef() =
        rho_.boundaryField()*(eu_.boundaryField() + K_.boundaryField());
}


void Foam::psiuCompressibleSystem::clear()
{
    rhoAdvection_.clear();
    rhoUAdvection_.clear();
    rhoEAdvection_.clear();
    rhoEuAdvection_.clear();
    rhoftAdvection_.clear();
    rhobAdvection_.clear();
    rhoSuAdvection_.clear();
    rhoXiAdvection_.clear();

    fluxScheme_->clear();
}


void Foam::psiuCompressibleSystem::encode()
{
    K_ = 0.5*magSqr(U_);
    rhoU_ = rho_*U_;
    rhoE_ = rho_*(e_ + K_);
    rhoEu_ = rho_*(eu_ + K_);
}


Foam::tmp<Foam::volScalarField>
Foam::psiuCompressibleSystem::speedOfSound() const
{
    return sqrt(thermo_->gamma()/thermo_->psi());
}


// ************************************************************************* //
