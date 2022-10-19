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

#include "XiCompressibleSystem.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(XiCompressibleSystem, 0);
    addToRunTimeSelectionTable
    (
        compressibleSystem,
        XiCompressibleSystem,
        singlePhase
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiCompressibleSystem::XiCompressibleSystem
(
    const fvMesh& mesh
)
:
    compressibleSystem(mesh),
    thermo_(psiuReactionThermo::New(mesh)),
    rho_
    (
        IOobject
        (
            "rho",
            mesh.time().timeName(),
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
            mesh.time().timeName(),
            mesh
        ),
        eu_*rho_
    ),
    composition_(thermo_->composition()),
    b_(composition_.Y("b")),
    Xi_
    (
        IOobject
        (
            "Xi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    unstrainedLaminarFlameSpeed(laminarFlameSpeed::New(thermo)),
    Su_
    (
        IOobject
        (
            "Su",
            runTime.timeName(),
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
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Xi*Su
    ),
    SuMin_(0.01*Su.average()),
    SuMax_(4.0*Su.average())
{
    thermo_->validate("XiCompressibleSystem", "ea");

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
        thermophysicalTransport_ =
            fluidThermophysicalTransportModel::New
            (
                turbulence_(),
                thermo_()
            );
    }

    fluxScheme_ = fluxScheme::NewSingle(mesh);
    encode();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::XiCompressibleSystem::~XiCompressibleSystem()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::XiCompressibleSystem::solve()
{
    volScalarField deltaRho(fvc::div(rhoPhi_));
    volVectorField deltaRhoU(fvc::div(rhoUPhi_) - g_*rho_);
    volScalarField deltaRhoE
    (
        fvc::div(rhoEPhi_)
      - (rhoU_ & g_)
    );
    volScalarField deltaRhoEu
    (
        fvc::div(fluxScheme_->energyFlux(rho_, U_, eu_, p_))
      - (rhoU_ & g_)
    );
    volScalarField deltab(fvc::div(rhoPhi_, b_));

    //- Store changed in mass, momentum and energy
    this->storeAndBlendDelta(deltaRho);
    this->storeAndBlendDelta(deltaRhoU);
    this->storeAndBlendDelta(deltaRhoE);
    this->storeAndBlendDelta(deltaRhoEu);
    this->storeAndBlendDelta(deltab);

    //- Store old values
    this->storeAndBlendOld(rho_);
    rho_.storePrevIter();

    this->storeAndBlendOld(rhoU_);
    this->storeAndBlendOld(rhoE_);
    this->storeAndBlendOld(rhoEu_);

    this->storeAndBlendOld(b_, false);

    dimensionedScalar dT = rho_.time().deltaT();

    rho_ -= dT*deltaRho;
    rho_.correctBoundaryConditions();

    vector solutionDs((vector(rho_.mesh().solutionD()) + vector::one)/2.0);
    rhoU_ -= cmptMultiply(dT*deltaRhoU, solutionDs);
    rhoE_ -= dT*deltaRhoE;
    rhoEu_ -= dT*deltaRhoEu;

    b_ = (b_*rho_.prevIter() - dT*deltab)/rho_;

    if (composition_.contains("ft"))
    {
        volScalarField& ft = composition_.Y("ft");

        volScalarField deltaft(fvc::div(rhoPhi_, ft_));
        this->storeAndBlendDelta(deltaft);

        this->storeAndBlendOld(ft_, false);

        ft = (ft*rho_.prevIter() - dT*deltaft)/rho_;
    }
}


void Foam::XiCompressibleSystem::postUpdate()
{
    this->decode();

   // Solve momentum diffusion
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho_, U_) - fvc::ddt(rho_, U_)
      + turbulence_->divDevTau(U_)
    );
    volScalarField dTDivSigmaDotU
    (
        rho_.mesh().time().deltaT()
       *fvc::div
        (
            fvc::dotInterpolate(rho_.mesh().Sf(), turbulence_->devTau())
          & fluxScheme_->Uf()
        )
    );
    rhoE_ += dTDivSigmaDotU;
    rhoEu_ += dTDivSigmaDotU;

    UEqn.solve();
    rhoU_ = rho_*U_;

    // Solve thermal energy diffusion
    e_ = rhoE_/rho_ - 0.5*magSqr(U_);
    eu_ = rhoEu_/rho_ - 0.5*magSqr(U_);
    Foam::solve
    (
        fvm::ddt(rho_, e_) - fvc::ddt(rho_, e_)
      - fvm::laplacian(thermophysicalTransport_->alphaEff(), e_)
    );
    Foam::solve
    (
        fvm::ddt(rho_, eu_) - fvc::ddt(rho_, eu_)
      - fvm::laplacian(thermophysicalTransport_->alphaEff(), eu_)
    );

    rhoE_ = rho_*(e_ + 0.5*magSqr(U_));
    rhoEu_ = rho_*(eu_ + 0.5*magSqr(U_));

    turbulence_->correct();

//     if(temperatureFix)
    {
        scalar Tulow = 250.0;
        volScalarField dummyTu(thermo_->Tu());
        dummyTu.min(Tulow);

        eu_ = max(eu_, thermo_->he(p_, dummyTu));
        rhoEu_ = rho_*(eu_ + 0.5*magSqr(U_));
    }

    thermo_->correct();
    p_.ref() = rho_/thermo_->psi();
    p_.correctBoundaryConditions();
    rho_.boundaryFieldRef() ==
        thermo_->psi().boundaryField()*p_.boundaryField();

    if (composition_.contains("ft"))
    {
        volScalarField& ft = composition_.Y("ft");
        fvScalarMatrix ftEqn
        (
            fvm::ddt(rho_, ft)
          - fvc::ddt(rho_, ft)
          - fvm::laplacian(thermophysicalTransport.alphaEff(), ft)
        );
        ftEqn.solve();
    }

    if (ign_.ignited())
    {
        // progress variable
        // ~~~~~~~~~~~~~~~~~
        volScalarField c("c", scalar(1) - b_);

        // Unburnt gas density
        // ~~~~~~~~~~~~~~~~~~~
        volScalarField rhou(thermo_->rhou());

        // Calculate flame normal etc.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        volVectorField n("n", fvc::grad(b_));

        volScalarField mgb(mag(n));

        dimensionedScalar dMgb =
            1.0e-3
           *(b_*c*mgb)().weightedAverage(mesh_.V())
           /((b_*c)().weightedAverage(mesh_.V()) + small)
          + dimensionedScalar(mgb.dimensions(), small);

        mgb += dMgb;

        surfaceVectorField SfHat(mesh_.Sf()/mesh_.magSf());
        surfaceVectorField nfVec(fvc::interpolate(n));
        nfVec += SfHat*(fvc::snGrad(b_) - (SfHat & nfVec));
        nfVec /= (mag(nfVec) + dMgb);
        surfaceScalarField nf((mesh_.Sf() & nfVec));
        n /= mgb;


        #include "StCorr.H"

        // Calculate turbulent flame speed flux
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        surfaceScalarField phiSt("phiSt", fvc::interpolate(rhou*StCorr*Su_*Xi_)*nf);

        scalar StCoNum = max
        (
            mesh_.surfaceInterpolation::deltaCoeffs()
           *mag(phiSt)/(fvc::interpolate(rho_)*mesh_.magSf())
        ).value()*mesh_.time().deltaTValue();

        Info<< "Max St-Courant Number = " << StCoNum << endl;

        // Create b equation
        // ~~~~~~~~~~~~~~~~~
        fvScalarMatrix bEqn
        (
            fvm::ddt(rho_, b_) - fvc::ddt(rho_, b_)
          + fvm::div(phiSt, b_)
          - fvm::Sp(fvc::div(phiSt), b)
          - fvm::laplacian(thermophysicalTransport.alphaEff(), b)
        );


        // Add ignition cell contribution to b-equation
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #include "ignite.H"


        // Solve for b
        // ~~~~~~~~~~~
        bEqn.solve();


        Info<< "min(b) = " << min(b_).value() << endl;


        // Calculate coefficients for Gulder's flame speed correlation
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        volScalarField up(uPrimeCoef*sqrt((2.0/3.0)*turbulence_->k()));
      // volScalarField up(sqrt(mag(diag(n * n) & diag(turbulence->r()))));

        volScalarField epsilon(pow(uPrimeCoef, 3)*turbulence_->epsilon());

        volScalarField tauEta(sqrt(thermo.muu()/(rhou*epsilon)));

        volScalarField Reta
        (
            up
          / (
                sqrt(epsilon*tauEta)
              + dimensionedScalar(up.dimensions(), 1e-8)
            )
        );

      // volScalarField l = 0.337*k*sqrt(k)/epsilon;
      // Reta *= max((l - dimensionedScalar(dimLength, 1.5e-3))/l, 0);

        // Calculate Xi flux
        // ~~~~~~~~~~~~~~~~~
        surfaceScalarField phiXi
        (
            phiSt
          - fvc::interpolate(fvc::laplacian(thermophysicalTransport.alphaEff(), b)/mgb)*nf
          + fvc::interpolate(rho)*fvc::interpolate(Su*(1.0/Xi - Xi))*nf
        );


        // Calculate mean and turbulent strain rates
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        volVectorField Ut(U + Su*Xi*n);
        volScalarField sigmat((n & n)*fvc::div(Ut) - (n & fvc::grad(Ut) & n));

        volScalarField sigmas
        (
            ((n & n)*fvc::div(U) - (n & fvc::grad(U) & n))/Xi
          + (
                (n & n)*fvc::div(Su*n)
              - (n & fvc::grad(Su*n) & n)
            )*(Xi + scalar(1))/(2*Xi)
        );


        // Calculate the unstrained laminar flame speed
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        volScalarField Su0(unstrainedLaminarFlameSpeed()());


        // Calculate the laminar flame speed in equilibrium with the applied strain
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        volScalarField SuInf(Su0*max(scalar(1) - sigmas/sigmaExt, scalar(0.01)));

        if (SuModel == "unstrained")
        {
            Su == Su0;
        }
        else if (SuModel == "equilibrium")
        {
            Su == SuInf;
        }
        else if (SuModel == "transport")
        {
            // Solve for the strained laminar flame speed
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            volScalarField Rc
            (
                (sigmas*SuInf*(Su0 - SuInf) + sqr(SuMin)*sigmaExt)
                /(sqr(Su0 - SuInf) + sqr(SuMin))
            );

            fvScalarMatrix SuEqn
            (
                fvm::ddt(rho, Su)
              + fvm::div(phi, Su)
              + fvm::div(phiXi, Su, "div(phiXi,Su)")
              - fvm::Sp(fvc::div(phiXi), Su)
              ==
              - fvm::SuSp(-rho*Rc*Su0/Su, Su)
              - fvm::SuSp(rho*(sigmas + Rc), Su)
            );
            SuEqn.solve();


            // Limit the maximum Su
            // ~~~~~~~~~~~~~~~~~~~~
            Su.min(SuMax);
            Su.max(SuMin);
        }
        else
        {
            FatalError
                << args.executable() << " : Unknown Su model " << SuModel
                << abort(FatalError);
        }


        // Calculate Xi according to the selected flame wrinkling model
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (XiModel == "fixed")
        {
            // Do nothing, Xi is fixed!
        }
        else if (XiModel == "algebraic")
        {
            // Simple algebraic model for Xi based on Gulders correlation
            // with a linear correction function to give a plausible profile for Xi
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            Xi == scalar(1) +
                (scalar(1) + (2*XiShapeCoef)*(scalar(0.5) - b))
               *XiCoef*sqrt(up/(Su + SuMin))*Reta;
        }
        else if (XiModel == "transport")
        {
            // Calculate Xi transport coefficients based on Gulders correlation
            // and DNS data for the rate of generation
            // with a linear correction function to give a plausible profile for Xi
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            volScalarField XiEqStar
            (
                scalar(1.001) + XiCoef*sqrt(up/(Su + SuMin))*Reta
            );

            volScalarField XiEq
            (
                scalar(1.001)
              + (
                    scalar(1)
                  + (2*XiShapeCoef)
                   *(scalar(0.5) - min(max(b, scalar(0)), scalar(1)))
                )*(XiEqStar - scalar(1.001))
            );

            volScalarField Gstar(0.28/tauEta);
            volScalarField R(Gstar*XiEqStar/(XiEqStar - scalar(1)));
            volScalarField G(R*(XiEq - scalar(1.001))/XiEq);

            // R *= (Gstar + 2*mag(dev(symm(fvc::grad(U)))))/Gstar;

            // Solve for the flame wrinkling
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            fvScalarMatrix XiEqn
            (
                fvm::ddt(rho, Xi)
              + fvm::div(phi, Xi)
              + fvm::div(phiXi, Xi, "div(phiXi,Xi)")
              - fvm::Sp(fvc::div(phiXi), Xi)
             ==
                rho*R
              - fvm::Sp(rho*(R - G), Xi)
              - fvm::Sp
                (
                    rho*max
                    (
                        sigmat - sigmas,
                        dimensionedScalar(sigmat.dimensions(), 0)
                    ),
                    Xi
                )
            );

            XiEqn.solve();

            // Correct boundedness of Xi
            // ~~~~~~~~~~~~~~~~~~~~~~~~~
            Xi.max(1.0);
            Info<< "max(Xi) = " << max(Xi).value() << endl;
            Info<< "max(XiEq) = " << max(XiEq).value() << endl;
        }
        else
        {
            FatalError
                << args.executable() << " : Unknown Xi model " << XiModel
                << abort(FatalError);
        }

        Info<< "Combustion progress = "
            << 100*(scalar(1) - b)().weightedAverage(mesh.V()).value() << "%"
            << endl;

        St = Xi*Su;
    }
}


void Foam::XiCompressibleSystem::update()
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


void Foam::XiCompressibleSystem::decode()
{
    U_.ref() = rhoU_()/rho_();
    U_.correctBoundaryConditions();

    rhoU_.boundaryFieldRef() = rho_.boundaryField()*U_.boundaryField();

    volScalarField E(rhoE_/rho_);
    e_.ref() = E() - 0.5*magSqr(U_());
    e_.correctBoundaryConditions();

    volScalarField Eu(rhoEu_/rho_);
    eu_.ref() = Eu() - 0.5*magSqr(U_());
    eu_.correctBoundaryConditions();

    rhoE_.boundaryFieldRef() =
        rho_.boundaryField()
       *(
            e_.boundaryField()
          + 0.5*magSqr(U_.boundaryField())
        );
    rhoEu_.boundaryFieldRef() =
        rho_.boundaryField()
       *(
            eu_.boundaryField()
          + 0.5*magSqr(U_.boundaryField())
        );

    thermo_->correct();
    p_.ref() = rho_/thermo_->psi();
    p_.correctBoundaryConditions();
    rho_.boundaryFieldRef() ==
        thermo_->psi().boundaryField()*p_.boundaryField();
}


void Foam::XiCompressibleSystem::encode()
{
    rhoU_ = rho_*U_;
    rhoE_ = rho_*(e_ + 0.5*magSqr(U_));
    rhoEu_ = rho_*(eu_ + 0.5*magSqr(U_));
}


Foam::tmp<Foam::volScalarField>
Foam::XiCompressibleSystem::speedOfSound() const
{
    return sqrt(thermo_->gamma()/thermo_->psi());
}


const Foam::volScalarField& Foam::XiCompressibleSystem::rho() const
{
    return rho_;
}


Foam::tmp<Foam::volScalarField> Foam::XiCompressibleSystem::rhou() const
{
    return thermo_->rhou();
}


Foam::tmp<Foam::volScalarField> Foam::XiCompressibleSystem::Cv() const
{
    return thermo_->Cv();
}


Foam::tmp<Foam::volScalarField> Foam::XiCompressibleSystem::mu() const
{
    return thermo_->mu();
}


Foam::tmp<Foam::scalarField>
Foam::XiCompressibleSystem::mu(const label patchi) const
{
    return thermo_->mu(patchi);
}


Foam::tmp<Foam::volScalarField> Foam::XiCompressibleSystem::muu() const
{
    return thermo_->muu();
}


Foam::tmp<Foam::volScalarField> Foam::XiCompressibleSystem::mub() const
{
    return thermo_->mub();
}


Foam::tmp<Foam::volScalarField> Foam::XiCompressibleSystem::nu() const
{
    return thermo_->nu();
}

Foam::tmp<Foam::scalarField>
Foam::XiCompressibleSystem::nu(const label patchi) const
{
    return thermo_->nu(patchi);
}

Foam::tmp<Foam::volScalarField>
Foam::XiCompressibleSystem::alpha() const
{
    return thermo_->alpha();
}

Foam::tmp<Foam::scalarField>
Foam::XiCompressibleSystem::alpha(const label patchi) const
{
    return thermo_->alpha(patchi);
}

Foam::tmp<Foam::volScalarField> Foam::XiCompressibleSystem::alphaEff
(
    const volScalarField& alphat
) const
{
    return thermo_->alphaEff(alphat);
}

Foam::tmp<Foam::scalarField> Foam::XiCompressibleSystem::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return thermo_->alphaEff(alphat, patchi);
}

Foam::tmp<Foam::volScalarField>
Foam::XiCompressibleSystem::alphahe() const
{
    return thermo_->alphahe();
}

Foam::tmp<Foam::scalarField>
Foam::XiCompressibleSystem::alphahe(const label patchi) const
{
    return thermo_->alphahe(patchi);
}

Foam::tmp<Foam::volScalarField> Foam::XiCompressibleSystem::kappa() const
{
    return thermo_->kappa();
}

Foam::tmp<Foam::scalarField>
Foam::XiCompressibleSystem::kappa(const label patchi) const
{
    return thermo_->kappa(patchi);
}

Foam::tmp<Foam::volScalarField> Foam::XiCompressibleSystem::kappaEff
(
    const volScalarField& alphat
) const
{
    return thermo_->kappaEff(alphat);
}

Foam::tmp<Foam::scalarField> Foam::XiCompressibleSystem::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return thermo_->kappaEff(alphat, patchi);
}
// ************************************************************************* //
