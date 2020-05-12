/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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

#include "coupledMultiphaseCompressibleSystem.H"
#include "fiveEqnCompressibleTurbulenceModel.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coupledMultiphaseCompressibleSystem, 0);
    addToRunTimeSelectionTable
    (
        phaseCompressibleSystem,
        coupledMultiphaseCompressibleSystem,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupledMultiphaseCompressibleSystem::coupledMultiphaseCompressibleSystem
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    multiphaseCompressibleSystem(mesh, dict),
    volumeFraction_
    (
        IOobject
        (
            "continuousVolumeFraction",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        0.0
    ),
    alphaRho_
    (
        IOobject
        (
            "alphaRho",
            mesh.time().timeName(),
            mesh
        ),
        volumeFraction_*rho_
    )
{
    forAll(alphas_, phasei)
    {
        volumeFraction_ += alphas_[phasei];
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coupledMultiphaseCompressibleSystem::~coupledMultiphaseCompressibleSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledMultiphaseCompressibleSystem::solve
(
    const label stepi,
    const scalarList& ai,
    const scalarList& bi
)
{
    PtrList<volScalarField> alphasOld(alphas_.size());
    PtrList<volScalarField> alphaRhosOld(alphas_.size());
    forAll(alphas_, phasei)
    {
        alphasOld.set
        (
            phasei, new volScalarField(alphas_[phasei])
        );
        alphaRhosOld.set
        (
            phasei, new volScalarField(alphaRhos_[phasei])
        );
    }
    if (oldIs_[stepi - 1] != -1)
    {
        rhoUOld_.set
        (
            oldIs_[stepi - 1],
            new volVectorField(rhoU_)
        );
        rhoEOld_.set
        (
            oldIs_[stepi - 1],
            new volScalarField(rhoE_)
        );
        forAll(alphas_, phasei)
        {
            alphasOld_[oldIs_[stepi - 1]].set
            (
                phasei,
                new volScalarField(alphas_[phasei])
            );
            alphaRhosOld_[oldIs_[stepi - 1]].set
            (
                phasei,
                new volScalarField(alphaRhos_[phasei])
            );
        }
    }

    volVectorField rhoUOld(ai[stepi - 1]*rhoU_);
    volScalarField rhoEOld(ai[stepi - 1]*rhoE_);
    forAll(alphas_, phasei)
    {
        alphasOld[phasei] *= ai[stepi - 1];
        alphaRhosOld[phasei] *= ai[stepi - 1];
    }
    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = oldIs_[i];
        if (fi != -1 && ai[fi] != 0)
        {
            rhoUOld += ai[fi]*rhoUOld_[fi];
            rhoEOld += ai[fi]*rhoEOld_[fi];
            forAll(alphas_, phasei)
            {
                alphasOld[phasei] += ai[fi]*alphasOld_[fi][phasei];
                alphaRhosOld[phasei] += ai[fi]*alphaRhosOld_[fi][phasei];
            }
        }
    }
    forAll(alphas_, phasei)
    {
        alphaRhos_[phasei].oldTime() = alphaRhosOld[phasei];
    }
    thermo_.solve(stepi, ai, bi);

    surfaceScalarField alphaf(fluxScheme_->interpolate(volumeFraction_, "alpha"));
    PtrList<volScalarField> deltaAlphas(alphas_.size());
    PtrList<volScalarField> deltaAlphaRhos(alphas_.size());
    volVectorField deltaRhoU
    (
        fvc::div(rhoUPhi_*alphaf)
      - p_*fvc::grad(alphaf)
      - g_*alphaRho_
    );
    volScalarField deltaRhoE
    (
        fvc::div(rhoEPhi_*alphaf)
      - ESource()
      - (rhoU_ & g_)
    );

    forAll(alphas_, phasei)
    {
        deltaAlphas.set
        (
            phasei,
            new volScalarField
            (
                fvc::div(alphaPhis_[phasei])
              - alphas_[phasei]*fvc::div(phi_)
            )
        );
        deltaAlphaRhos.set
        (
            phasei, new volScalarField(fvc::div(alphaRhoPhis_[phasei]))
        );
    }
    if (deltaIs_[stepi - 1] != -1)
    {
        deltaRhoU_.set
        (
            deltaIs_[stepi - 1],
            new volVectorField(deltaRhoU)
        );
        deltaRhoE_.set
        (
            deltaIs_[stepi - 1],
            new volScalarField(deltaRhoE)
        );
        forAll(alphas_, phasei)
        {
            deltaAlphas_[deltaIs_[stepi - 1]].set
            (
                phasei,
                new volScalarField(deltaAlphas[phasei])
            );
            deltaAlphaRhos_[deltaIs_[stepi - 1]].set
            (
                phasei,
                new volScalarField(deltaAlphaRhos[phasei])
            );
        }
    }
    deltaRhoU *= bi[stepi - 1];
    deltaRhoE *= bi[stepi - 1];
    forAll(alphas_, phasei)
    {
        deltaAlphas[phasei] *= bi[stepi - 1];
        deltaAlphaRhos[phasei] *= bi[stepi - 1];
    }

    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = deltaIs_[i];
        if (fi != -1 && bi[fi] != 0)
        {
            deltaRhoU += bi[fi]*deltaRhoU_[fi];
            deltaRhoE += bi[fi]*deltaRhoE_[fi];
            forAll(alphas_, phasei)
            {
                deltaAlphas[phasei] += bi[fi]*deltaAlphas_[fi][phasei];
                deltaAlphaRhos[phasei] += bi[fi]*deltaAlphaRhos_[fi][phasei];
            }
        }
    }

    dimensionedScalar dT = rho_.time().deltaT();
    vector solutionDs((vector(rho_.mesh().solutionD()) + vector::one)/2.0);

    rhoU_ = cmptMultiply(rhoUOld - dT*deltaRhoU, solutionDs);
    rhoE_ = rhoEOld - dT*deltaRhoE;
//     if (radiation_->type() != "none")
//     {
//         calcAlphaAndRho();
//         e() = rhoE_/alphaRho_ - 0.5*magSqr(U_);
//         e().correctBoundaryConditions();
//         rhoE_ = radiation_->calcRhoE(f*dT, rhoE_, alphaRho_, e(), Cv());
//     }

    forAll(alphas_, phasei)
    {
        alphas_[phasei] = alphasOld[phasei] - dT*deltaAlphas[phasei];
        alphas_[phasei].correctBoundaryConditions();

        alphaRhos_[phasei].oldTime() = alphaRhosOld[phasei];
        alphaRhos_[phasei] = alphaRhosOld[phasei] - dT*deltaAlphaRhos[phasei];
        alphaRhos_[phasei].correctBoundaryConditions();
    }


//     if (stepi == oldIs_.size())
//     {
//         radiation_->correct();
//     }

    if
    (
        stepi == oldIs_.size()
     && (
            turbulence_.valid()
         || dragSource_.valid()
         || extESource_.valid()
        )
    )
    {
        calcAlphaAndRho();
        U_ = rhoU_/alphaRho_;
        U_.correctBoundaryConditions();

        e() = rhoE_/alphaRho_ - 0.5*magSqr(U_);
        e().correctBoundaryConditions();

        fvVectorMatrix UEqn
        (
            fvm::ddt(alphaRho_, U_) - fvc::ddt(alphaRho_, U_)
        );
        fvScalarMatrix eEqn
        (
            fvm::ddt(alphaRho_, e()) - fvc::ddt(alphaRho_, e())
        );

        if (turbulence_.valid())
        {
            volScalarField muEff
            (
                "muEff",
                volumeFraction_*turbulence_->muEff()
            );
            UEqn -=
                fvm::laplacian(muEff, U_)
              + fvc::div(muEff*dev2(Foam::T(fvc::grad(U_))));

            eEqn -=
                fvm::laplacian
                (
                    volumeFraction_*turbulence_->alphaEff(),
                    e()
                );
        }

        if (dragSource_.valid())
        {
            UEqn -= dragSource_();
        }
        if (extESource_.valid())
        {
            eEqn -= extESource_();
        }
        UEqn.solve();
        eEqn.solve();

        rhoU_ = alphaRho_*U_;
        rhoE_ = alphaRho_*(e() + 0.5*magSqr(U_)); // Includes change to total energy from viscous term in momentum equation

        if (turbulence_.valid())
        {
            turbulence_->correct();
        }
    }

    decode();
}

void Foam::coupledMultiphaseCompressibleSystem::calcAlphaAndRho()
{
    alphaRho_ = dimensionedScalar("0", dimDensity, 0.0);
    volScalarField sumAlpha
    (
        IOobject
        (
            "sumAlpha",
            rho_.time().timeName(),
            rho_.mesh()
        ),
        rho_.mesh(),
        0.0
    );
    for (label phasei = 0; phasei < alphas_.size() - 1; phasei++)
    {
        alphas_[phasei].max(0);
        alphas_[phasei].min(1);
        alphas_[phasei].correctBoundaryConditions();
        sumAlpha += alphas_[phasei];

        alphaRhos_[phasei].max(0);
        rhos_[phasei] =
            alphaRhos_[phasei]
           /max(alphas_[phasei], thermo_.thermo(phasei).residualAlpha());
        rhos_[phasei].correctBoundaryConditions();

        alphaRhos_[phasei] = alphas_[phasei]*rhos_[phasei];
        alphaRho_ += alphaRhos_[phasei];
    }
    label phasei = alphas_.size() - 1;
    alphas_[phasei] = volumeFraction_ - sumAlpha;
    alphas_[phasei].max(0.0);
    alphas_[phasei].min(1.0);

    alphaRhos_[phasei].max(0.0);
    rhos_[phasei] = alphaRhos_[phasei]
           /max(alphas_[phasei], thermo_.thermo(phasei).residualAlpha());
    rhos_[phasei].correctBoundaryConditions();
    alphaRhos_[phasei] = alphas_[phasei]*rhos_[phasei];
    alphaRho_ += alphaRhos_[phasei];
    rho_ = alphaRho_/volumeFraction_;

}


void Foam::coupledMultiphaseCompressibleSystem::decode()
{
    calcAlphaAndRho();

    U_.ref() = rhoU_()/alphaRho_();
    U_.correctBoundaryConditions();

    rhoU_.boundaryFieldRef() = alphaRho_.boundaryField()*U_.boundaryField();

    volScalarField E(rhoE_/alphaRho_);
    e_.ref() = E() - 0.5*magSqr(U_());

    //--- Hard limit, e
    if(min(e_).value() < 0)
    {
        WarningInFunction<< "Limiting e, min(e) = " << min(e_).value() << endl;
        e_.max(small);
        rhoE_.ref() = alphaRho_()*(e_() + 0.5*magSqr(U_()));
    }
    e_.correctBoundaryConditions();

    rhoE_.boundaryFieldRef() =
        alphaRho_.boundaryField()
       *(
            e_.boundaryField()
          + 0.5*magSqr(U_.boundaryField())
        );

    thermo_.correct();
    T_ /= Foam::max(volumeFraction_, 1e-10);
}


void Foam::coupledMultiphaseCompressibleSystem::encode()
{
    alphaRho_ = dimensionedScalar("0", dimDensity, 0.0);
    forAll(alphas_, phasei)
    {
        alphaRhos_[phasei] = alphas_[phasei]*rhos_[phasei];
        alphaRho_ += alphaRhos_[phasei];
    }
    rho_ = alphaRho_/volumeFraction_;
    rhoU_ = alphaRho_*U_;
    rhoE_ = alphaRho_*(e_ + 0.5*magSqr(U_));
}



Foam::tmp<Foam::volScalarField> Foam::coupledMultiphaseCompressibleSystem::Cv() const
{
    return thermo_.Cv()/Foam::max(volumeFraction_, 1e-10);
}


Foam::tmp<Foam::volScalarField> Foam::coupledMultiphaseCompressibleSystem::mu() const
{
    return thermo_.mu()/Foam::max(volumeFraction_, 1e-10);
}


Foam::tmp<Foam::scalarField>
Foam::coupledMultiphaseCompressibleSystem::mu(const label patchi) const
{
    return
    thermo_.mu(patchi)
       /Foam::max(volumeFraction_.boundaryField()[patchi], 1e-10);
}

Foam::tmp<Foam::volScalarField> Foam::coupledMultiphaseCompressibleSystem::nu() const
{
    return thermo_.nu()/Foam::max(volumeFraction_, 1e-10);
}

Foam::tmp<Foam::scalarField>
Foam::coupledMultiphaseCompressibleSystem::nu(const label patchi) const
{
    return
        thermo_.nu(patchi)
       /Foam::max(volumeFraction_.boundaryField()[patchi], 1e-10);
}

Foam::tmp<Foam::volScalarField>
Foam::coupledMultiphaseCompressibleSystem::alpha() const
{
    return thermo_.alpha()/Foam::max(volumeFraction_, 1e-10);
}

Foam::tmp<Foam::scalarField>
Foam::coupledMultiphaseCompressibleSystem::alpha(const label patchi) const
{
    return
        thermo_.alpha(patchi)
       /Foam::max(volumeFraction_.boundaryField()[patchi], 1e-10);
}

Foam::tmp<Foam::volScalarField> Foam::coupledMultiphaseCompressibleSystem::alphaEff
(
    const volScalarField& alphat
) const
{
    return thermo_.alphaEff(alphat)/Foam::max(volumeFraction_, 1e-10);
}

Foam::tmp<Foam::scalarField> Foam::coupledMultiphaseCompressibleSystem::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        thermo_.alphaEff(alphat, patchi)
       /Foam::max(volumeFraction_.boundaryField()[patchi], 1e-10);
}

Foam::tmp<Foam::volScalarField>
Foam::coupledMultiphaseCompressibleSystem::alphahe() const
{
    return thermo_.alphahe()/Foam::max(volumeFraction_, 1e-10);
}

Foam::tmp<Foam::scalarField>
Foam::coupledMultiphaseCompressibleSystem::alphahe(const label patchi) const
{
    return
        thermo_.alphahe(patchi)
       /Foam::max(volumeFraction_.boundaryField()[patchi], 1e-10);
}

Foam::tmp<Foam::volScalarField> Foam::coupledMultiphaseCompressibleSystem::kappa() const
{
    return thermo_.kappa()/Foam::max(volumeFraction_, 1e-10);
}

Foam::tmp<Foam::scalarField>
Foam::coupledMultiphaseCompressibleSystem::kappa(const label patchi) const
{
    return
        thermo_.kappa(patchi)
       /Foam::max(volumeFraction_.boundaryField()[patchi], 1e-10);
}

Foam::tmp<Foam::volScalarField> Foam::coupledMultiphaseCompressibleSystem::kappaEff
(
    const volScalarField& alphat
) const
{
    return
        thermo_.kappaEff(alphat)
       /Foam::max(volumeFraction_, 1e-10);
}

Foam::tmp<Foam::scalarField> Foam::coupledMultiphaseCompressibleSystem::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        thermo_.kappaEff(alphat, patchi)
       /Foam::max(volumeFraction_.boundaryField()[patchi], 1e-10);
}

// ************************************************************************* //
