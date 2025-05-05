/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(psiuCompressibleSystem, 0);
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
            IOobject::READ_IF_PRESENT,
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

    if (min(thermo_->mu()).value() > small)
    {
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
    }

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

    //- Store changed in mass, momentum and energy
    this->storeAndBlendDelta(deltaRho);
    this->storeAndBlendDelta(deltaRhoU);
    this->storeAndBlendDelta(deltaRhoE);
    this->storeAndBlendDelta(deltaRhoEu);

    //- Store old values
    this->storeAndBlendOld(rho_);
    const volScalarField rho0(rho_);

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

    const volScalarField twoMinusf(2.0 - rho_/rho0);
    if (thermo_->containsSpecie("ft"))
    {
        volScalarField& ft = thermo_->Y("ft");
        volScalarField deltaRhoFt
        (
            fvc::div(fluxScheme_->flux(ft, rhoPhi_))
          - fvc::laplacian(thermophysicalTransport_->DEff(ft), ft)
        );
        this->storeAndBlendDelta(deltaRhoFt);
        this->storeAndBlendOld(ft);

        ft = ft*twoMinusf - dT*deltaRhoFt/rho0;
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
        const volScalarField& rho = rho_;

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
           *mag(phiSt)/(fvc::interpolate(rho)*mesh.magSf())
        ).value()*rho.time().deltaTValue();

        Info<< "Max St-Courant Number = " << StCoNum << endl;

        tmp<surfaceScalarField> tbf(fluxScheme_->interpolate(b_, "b"));
        const surfaceScalarField& bf = tbf();

        volScalarField deltaRhoB
        (
            "deltaRhoB",
            fvc::div(rhoPhi_*bf)
          + fvc::div(phiSt*bf)
          - b_*fvc::div(phiSt)
          - fvc::laplacian(thermophysicalTransport_->DEff(b_), b_)
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

                    deltaRhoB[ignCell] -=
                        (
                            ignSite.strength()*rhou[ignCell]/ignSite.duration()
                        )/(b_[ignCell] + 0.001);
                }
            }
        }

        const volScalarField bOld(b_);

        this->storeAndBlendDelta(deltaRhoB);
        this->storeAndBlendOld(b_);
        b_ = b_*twoMinusf - dT*deltaRhoB/rho0;
        b_.maxMin(0.0, 1.0);
        b_.correctBoundaryConditions();


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
              + fvc::interpolate(rho)*fvc::interpolate(Su_*(1.0/Xi_ - Xi_))*nf;
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
          + fvc::interpolate(rho)*fvc::interpolate(Su_*(1.0/Xi_ - Xi_))*nf
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

            volScalarField deltaRhoSu
            (
                "deltaRhoSu",
                fvc::div((phi_ + phiXi)*fluxScheme_->interpolate(Su_, "Su"))
              - Su_*(fvc::div(phiXi) + rho*(Rc*Su0 - (sigmas + Rc)))
            );
            this->storeAndBlendDelta(deltaRhoSu);
            this->storeAndBlendOld(Su_);
            Su_ = Su_*twoMinusf - dT*deltaRhoSu/rho0;
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
                  + (1.0 + (2*XiShapeCoef_)*(0.5 - bOld))
                   *XiCoef_*sqrt(up/(SuOld + SuMin_))*Reta;
            }
            else
            {
                volVectorField Ut(U_ + SuOld*Xi_*n);
                volScalarField sigmat(nDotn*fvc::div(Ut) - (n & fvc::grad(Ut) & n));

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
                 + (1.0 + (2*XiShapeCoef_)*(0.5 - bOld))*(XiEqStar - 1.001)
                );

                volScalarField Gstar(0.28/tauEta);
                volScalarField R(Gstar*XiEqStar/(XiEqStar - 1.0));
                volScalarField G(R*(XiEq - 1.001)/XiEq);

                // R *= (Gstar + 2*mag(dev(symm(fvc::grad(U)))))/Gstar;

                // Solve for the flame wrinkling
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                tmp<surfaceScalarField> tXif(fluxScheme_->interpolate(Xi_, "Xi"));
                const surfaceScalarField& Xif = tXif();

                volScalarField deltaRhoXi
                (
                    "deltaRhoXi",
                    fvc::div((phi_ + phiXi)*Xif)
                  - fvc::div(phiXi)*Xi_
                  + rho
                   *(
                        (Xi_ - 1.0)*R - G
                      + max
                        (
                            sigmat - tsigmas,
                            dimensionedScalar(sigmat.dimensions(), 0)
                        )*Xi_
                    )
                );
                this->storeAndBlendDelta(deltaRhoXi);
                this->storeAndBlendOld(Xi_);
                Xi_ = Xi_*twoMinusf - dT*deltaRhoXi/rho0;
                Xi_.max(1.0);
                Xi_.correctBoundaryConditions();

                Info<< "max(Xi) = " << max(Xi_).value() << endl;
                Info<< "max(XiEq) = " << max(XiEq).value() << endl;
            }
        }
        else
        {
            FatalError
                << "Unknown Xi model " << SuModel_
                << abort(FatalError);
        }
        St_ = Xi_*Su_;
    }
}


void Foam::psiuCompressibleSystem::postUpdate()
{
    if (!turbulence_.valid())
    {
        return;
    }

    this->decode();
    turbulence_->correct();
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

    rhoU_.boundaryFieldRef() = rho_.boundaryField()*U_.boundaryField();

    e_.internalFieldRef() = rhoE_()/rho_() - 0.5*magSqr(U_());
    e_.correctBoundaryConditions();
    rhoE_.boundaryFieldRef() =
        rho_.boundaryField()
       *(
            e_.boundaryField()
          + 0.5*magSqr(U_.boundaryField())
        );

    eu_.internalFieldRef() = rhoEu_()/rho_() - 0.5*magSqr(U_());
    forAll(b_, i)
    {
        if (b_[i] < 1e-6)
        {
            eu_[i] = e_[i];
        }
    }
    eu_.correctBoundaryConditions();
    rhoEu_.boundaryFieldRef() =
        rho_.boundaryField()
       *(
            eu_.boundaryField()
          + 0.5*magSqr(U_.boundaryField())
        );

    thermo_->correct();
    p_.internalFieldRef() = rho_/thermo_->psi();
    p_.correctBoundaryConditions();
    rho_.boundaryFieldRef() ==
        thermo_->psi().boundaryField()*p_.boundaryField();
}


void Foam::psiuCompressibleSystem::encode()
{
    rhoU_ = rho_*U_;
    rhoE_ = rho_*(e_ + 0.5*magSqr(U_));
    rhoEu_ = rho_*(eu_ + 0.5*magSqr(U_));
}


Foam::tmp<Foam::volScalarField>
Foam::psiuCompressibleSystem::speedOfSound() const
{
    return sqrt(thermo_->gamma()/thermo_->psi());
}


// ************************************************************************* //
