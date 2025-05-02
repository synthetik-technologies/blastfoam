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
    const fvMesh& mesh
)
:
    compressibleSystem(mesh),
    thermo_(psiuReactionThermo::New(mesh)),
    composition_(thermo_->composition()),
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

    combustionProperties
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

    b(composition_.Y("b")),
    Xi
    (
        IOobject
        (
            "Xi",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    Su
    (
        IOobject
        (
            "Su",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    St
    (
        IOobject
        (
            "St",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Xi*Su
    ),

    SuMin(0.01*Su.average()),
    SuMax(4.0*Su.average()),

    unstrainedLaminarFlameSpeed_(laminarFlameSpeed::New(thermo_())),

    SuModel(combustionProperties.lookup("SuModel")),
    XiModel(combustionProperties.lookup("XiModel")),

    sigmaExt(combustionProperties.lookup("sigmaExt")),
    XiCoef(combustionProperties.lookup("XiCoef")),
    XiShapeCoef(combustionProperties.lookup("XiShapeCoef")),
    uPrimeCoef(combustionProperties.lookup("uPrimeCoef")),

    ign(combustionProperties, mesh.time(), mesh)
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
        thermophysicalTransport_ =
            fluidThermophysicalTransportModel::New
            (
                turbulence_(),
                thermo_()
            );
    }

    fluxScheme_ = fluxScheme::NewSingle(phi_);
    encode();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::psiuCompressibleSystem::~psiuCompressibleSystem()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::psiuCompressibleSystem::solve()
{
    tmp<surfaceScalarField> talphaf
    (
        fvc::interpolate(thermophysicalTransport_->alphaEff())
    );
    const surfaceScalarField& alphaf = talphaf();

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
      - fvc::div(devTau)
    );
    volScalarField deltaRhoE
    (
        fvc::div(rhoEPhi_)
      - (rhoU_ & g_)
      - fvc::laplacian(alphaf, e_)
      - divSigmaDotU
    );
    volScalarField deltaRhoEu
    (
        fvc::div(fluxScheme_->energyFlux(rho_, U_, eu_, p_))
      - (rhoU_ & g_)
      - fvc::laplacian(alphaf, eu_)
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
    if (thermo_->composition().contains("ft"))
    {
        volScalarField& ft = thermo_->composition().Y("ft");
        volScalarField deltaRhoFt
        (
            fvc::div(fluxScheme_->flux(ft, rhoPhi_))
          - fvc::laplacian(alphaf, ft)
        );
        this->storeAndBlendDelta(deltaRhoFt);
        this->storeAndBlendOld(ft);

        ft = ft*twoMinusf - dT*deltaRhoFt/rho0;
        ft.max(0.0);
        ft.correctBoundaryConditions();
    }

    if (ign.ignited())
    {
        const fvMesh& mesh = this->mesh();

        // progress variable
        // ~~~~~~~~~~~~~~~~~
        volScalarField c("c", scalar(1.0) - b);

        // Unburnt gas density
        // ~~~~~~~~~~~~~~~~~~~
        volScalarField rhou(thermo_->rhou());
        const volScalarField& rho = rho_;

        // Calculate flame normal etc.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        volVectorField n("n", fvc::grad(b));

        volScalarField mgb(mag(n));

        dimensionedScalar dMgb =
            1.0e-3
           *(b*c*mgb)().weightedAverage(mesh.V())
           /((b*c)().weightedAverage(mesh.V()) + small)
          + dimensionedScalar(mgb.dimensions(), small);
        mgb += dMgb;

        surfaceVectorField SfHat(mesh.Sf()/mesh.magSf());
        surfaceVectorField nfVec(fvc::interpolate(n));
        nfVec += SfHat*(fvc::snGrad(b) - (SfHat & nfVec));
        nfVec /= (mag(nfVec) + dMgb);
        surfaceScalarField nf((mesh.Sf() & nfVec));
        n /= mgb;

        const volScalarField nDotn(n & n);



        #include "StCorr.H"

        // Calculate turbulent flame speed flux
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        surfaceScalarField phiSt("phiSt", fvc::interpolate(rhou*StCorr*Su*Xi)*nf);

        scalar StCoNum = max
        (
            mesh.surfaceInterpolation::deltaCoeffs()
           *mag(phiSt)/(fvc::interpolate(rho)*mesh.magSf())
        ).value()*rho.time().deltaTValue();

        Info<< "Max St-Courant Number = " << StCoNum << endl;

        tmp<surfaceScalarField> tbf(fluxScheme_->interpolate(b, "b"));
        const surfaceScalarField& bf = tbf();

        volScalarField deltaRhoB
        (
            "deltaRhoB",
            fvc::div(rhoPhi_*bf)
          + fvc::div(phiSt*bf)
          - b*fvc::div(phiSt)
          - fvc::laplacian(alphaf, b)
        );
        forAll(ign.sites(), i)
        {
            const ignitionSite& ignSite = ign.sites()[i];

            if (ignSite.igniting())
            {
                forAll(ignSite.cells(), icelli)
                {
                    label ignCell = ignSite.cells()[icelli];
                    DebugInfo<< "Igniting cell " << ignCell;

                    DebugInfo<< " state :"
                        << ' ' << b[ignCell]
                        << ' ' << Xi[ignCell]
                        << ' ' << Su[ignCell]
                        << ' ' << mgb[ignCell]
                        << endl;

                    deltaRhoB[ignCell] -=
                        (
                            ignSite.strength()*rhou[ignCell]/ignSite.duration()
                        )/(b[ignCell] + 0.001);
                }
            }
        }

        const volScalarField bOld(b);

        this->storeAndBlendDelta(deltaRhoB);
        this->storeAndBlendOld(b);
        b = b*twoMinusf - dT*deltaRhoB/rho0;
        b.maxMin(0.0, 1.0);
        b.correctBoundaryConditions();


        // Calculate Xi flux
        // ~~~~~~~~~~~~~~~~~
        tmp<surfaceScalarField> tphiXi;
        if (SuModel == "transport" || XiModel == "transport")
        {
            tphiXi =
                phiSt
              - fvc::interpolate(fvc::laplacian(alphaf, bOld)/mgb)*nf
              + fvc::interpolate(rho)*fvc::interpolate(Su*(1.0/Xi - Xi))*nf;
        }

        tmp<volScalarField> tsigmas;
        if
        (
            SuModel == "equilibrium"
         || SuModel == "transport"
         || XiModel == "transport"
        )
        (
            tsigmas =
                (nDotn*fvc::div(phi_) - (n & fvc::grad(U_) & n))/Xi
              + (
                    nDotn*fvc::div(Su*n)
                  - (n & fvc::grad(Su*n) & n)
                )*(Xi + scalar(1))/(2*Xi)
        );

        // Calculate Xi flux
        // ~~~~~~~~~~~~~~~~~
        surfaceScalarField phiXi
        (
            phiSt
          - fvc::interpolate(fvc::laplacian(alphaf, bOld)/mgb)*nf
          + fvc::interpolate(rho)*fvc::interpolate(Su*(1.0/Xi - Xi))*nf
        );

        const volScalarField SuOld(Su);

        volScalarField Su0(unstrainedLaminarFlameSpeed_()());
        if (SuModel == "unstrained")
        {
            Su == Su0;
        }
        else if (SuModel == "equilibrium")
        {
            Su == Su0*max(scalar(1) - tsigmas/sigmaExt, scalar(0.01));
        }
        else if (SuModel == "transport")
        {
            const volScalarField& sigmas = tsigmas();
            volScalarField SuInf(Su0*max(scalar(1) - sigmas/sigmaExt, scalar(0.01)));

            // Solve for the strained laminar flame speed
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            volScalarField Rc
            (
                (sigmas*SuInf*(Su0 - SuInf) + sqr(SuMin)*sigmaExt)
                /(sqr(Su0 - SuInf) + sqr(SuMin))
            );

            volScalarField deltaRhoSu
            (
                "deltaRhoSu",
                fvc::div((phi_ + phiXi)*fluxScheme_->interpolate(Su, "Su"))
              - Su*(fvc::div(phiXi) + rho*(Rc*Su0 - (sigmas + Rc)))
            );
            this->storeAndBlendDelta(deltaRhoSu);
            this->storeAndBlendOld(Su);
            Su = Su*twoMinusf - dT*deltaRhoSu/rho0;
            Su.maxMin(SuMin, SuMin);
            Su.correctBoundaryConditions();
        }
        else
        {
            FatalError
                << "Unknown Su model " << SuModel
                << abort(FatalError);
        }


        // Calculate Xi according to the selected flame wrinkling model
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (XiModel == "fixed")
        {
            // Do nothing, Xi is fixed!
        }
        else if (XiModel == "algebraic" || XiModel == "transport")
        {
            volScalarField epsilon(pow(uPrimeCoef, 3)*turbulence().epsilon());
            volScalarField tauEta(sqrt(thermo_->muu()/(rhou*epsilon)));

            volScalarField up(uPrimeCoef*sqrt((2.0/3.0)*turbulence().k()));
//             volScalarField up(sqrt(mag(diag(n * n) & diag(turbulence->r()))));

            tmp<volScalarField> Reta
            (
                up/(sqrt(epsilon*tauEta)+ dimensionedScalar(up.dimensions(), 1e-8))
            );
//             volScalarField l = 0.337*k*sqrt(k)/epsilon;
//             Reta *= max((l - dimensionedScalar(dimLength, 1.5e-3))/l, 0);

            if (XiModel == "algebraic")
            {
                // Simple algebraic model for Xi based on Gulders correlation
                // with a linear correction function to give a plausible profile
                // for Xi
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                Xi ==
                    1.0
                  + (1.0 + (2*XiShapeCoef)*(0.5 - bOld))
                   *XiCoef*sqrt(up/(SuOld + SuMin))*Reta;
            }
            else
            {
                volVectorField Ut(U_ + SuOld*Xi*n);
                volScalarField sigmat(nDotn*fvc::div(Ut) - (n & fvc::grad(Ut) & n));

                // Calculate Xi transport coefficients based on Gulders correlation
                // and DNS data for the rate of generation
                // with a linear correction function to give a plausible profile for Xi
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                volScalarField XiEqStar
                (
                    1.001 + XiCoef*sqrt(up/(SuOld + SuMin))*Reta
                );
                volScalarField XiEq
                (
                    1.001
                 + (1.0 + (2*XiShapeCoef)*(0.5 - bOld))*(XiEqStar - 1.001)
                );

                volScalarField Gstar(0.28/tauEta);
                volScalarField R(Gstar*XiEqStar/(XiEqStar - 1.0));
                volScalarField G(R*(XiEq - 1.001)/XiEq);

                // R *= (Gstar + 2*mag(dev(symm(fvc::grad(U)))))/Gstar;

                // Solve for the flame wrinkling
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                tmp<surfaceScalarField> tXif(fluxScheme_->interpolate(Xi, "Xi"));
                const surfaceScalarField& Xif = tXif();

                volScalarField deltaRhoXi
                (
                    "deltaRhoXi",
                    fvc::div((phi_ + phiXi)*Xif)
                  - fvc::div(phiXi)*Xi
                  + rho
                   *(
                        (Xi - 1.0)*R - G
                      + max
                        (
                            sigmat - tsigmas,
                            dimensionedScalar(sigmat.dimensions(), 0)
                        )*Xi
                    )
                );
                this->storeAndBlendDelta(deltaRhoXi);
                this->storeAndBlendOld(Xi);
                Xi = Xi*twoMinusf - dT*deltaRhoXi/rho0;
                Xi.max(1.0);
                Xi.correctBoundaryConditions();

                Info<< "max(Xi) = " << max(Xi).value() << endl;
                Info<< "max(XiEq) = " << max(XiEq).value() << endl;
            }
        }
        else
        {
            FatalError
                << "Unknown Xi model " << SuModel
                << abort(FatalError);
        }
        St = Xi*Su;
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
    U_.ref() = rhoU_()/rho_();
    U_.correctBoundaryConditions();

    rhoU_.boundaryFieldRef() = rho_.boundaryField()*U_.boundaryField();

    e_.ref() = rhoE_()/rho_() - 0.5*magSqr(U_());
    e_.correctBoundaryConditions();
    rhoE_.boundaryFieldRef() =
        rho_.boundaryField()
       *(
            e_.boundaryField()
          + 0.5*magSqr(U_.boundaryField())
        );

    eu_.ref() = rhoEu_()/rho_() - 0.5*magSqr(U_());
    forAll(b, i)
    {
        if (b[i] < 1e-6)
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
    p_.ref() = rho_/thermo_->psi();
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


Foam::tmp<Foam::volScalarField> Foam::psiuCompressibleSystem::rhou() const
{
    return thermo_->rhou();
}


Foam::tmp<Foam::volScalarField> Foam::psiuCompressibleSystem::Cv() const
{
    return thermo_->Cv();
}


Foam::tmp<Foam::volScalarField> Foam::psiuCompressibleSystem::mu() const
{
    return thermo_->mu();
}


Foam::tmp<Foam::scalarField>
Foam::psiuCompressibleSystem::mu(const label patchi) const
{
    return thermo_->mu(patchi);
}


Foam::tmp<Foam::volScalarField> Foam::psiuCompressibleSystem::muu() const
{
    return thermo_->muu();
}


Foam::tmp<Foam::volScalarField> Foam::psiuCompressibleSystem::mub() const
{
    return thermo_->mub();
}


Foam::tmp<Foam::volScalarField> Foam::psiuCompressibleSystem::nu() const
{
    return thermo_->nu();
}

Foam::tmp<Foam::scalarField>
Foam::psiuCompressibleSystem::nu(const label patchi) const
{
    return thermo_->nu(patchi);
}

Foam::tmp<Foam::volScalarField>
Foam::psiuCompressibleSystem::alpha() const
{
    return thermo_->alpha();
}

Foam::tmp<Foam::scalarField>
Foam::psiuCompressibleSystem::alpha(const label patchi) const
{
    return thermo_->alpha(patchi);
}

Foam::tmp<Foam::volScalarField> Foam::psiuCompressibleSystem::alphaEff
(
    const volScalarField& alphat
) const
{
    return thermo_->alphaEff(alphat);
}

Foam::tmp<Foam::scalarField> Foam::psiuCompressibleSystem::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return thermo_->alphaEff(alphat, patchi);
}

Foam::tmp<Foam::volScalarField>
Foam::psiuCompressibleSystem::alphahe() const
{
    return thermo_->alphahe();
}

Foam::tmp<Foam::scalarField>
Foam::psiuCompressibleSystem::alphahe(const label patchi) const
{
    return thermo_->alphahe(patchi);
}

Foam::tmp<Foam::volScalarField> Foam::psiuCompressibleSystem::kappa() const
{
    return thermo_->kappa();
}

Foam::tmp<Foam::scalarField>
Foam::psiuCompressibleSystem::kappa(const label patchi) const
{
    return thermo_->kappa(patchi);
}

Foam::tmp<Foam::volScalarField> Foam::psiuCompressibleSystem::kappaEff
(
    const volScalarField& alphat
) const
{
    return thermo_->kappaEff(alphat);
}

Foam::tmp<Foam::scalarField> Foam::psiuCompressibleSystem::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return thermo_->kappaEff(alphat, patchi);
}
// ************************************************************************* //
