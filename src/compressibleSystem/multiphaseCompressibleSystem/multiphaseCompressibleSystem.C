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

#include "multiphaseCompressibleSystem.H"
#include "addToRunTimeSelectionTable.H"
#include "SortableList.H"
#include "MULES.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiphaseCompressibleSystem, 0);
    addToRunTimeSelectionTable
    (
        compressibleSystem,
        multiphaseCompressibleSystem,
        multiphase
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseCompressibleSystem::multiphaseCompressibleSystem
(
    const fvMesh& mesh,
    const bool initialize
)
:
    compressibleBlastSystem(mesh, multiphaseFluidBlastThermo::typeName),
    thermo_(dynamicCast<multiphaseFluidBlastThermo>(thermoPtr_())),
    alphas_(thermo_.volumeFractions()),
    rhos_(thermo_.rhos()),
    alphaRhos_(alphas_.size()),
    alphaPhis_(alphas_.size()),
    alphaRhoPhis_(alphas_.size()),
    transportPhaseDensity_(this->lookupOrDefault("transportPhaseDensity", false)),
    sharpen_(alphas_.size(), false)

{
    this->fluxScheme_ = fluxScheme::NewMulti(phi_);

    bool allRead = rhoU_.headerOk() && rhoE_.headerOk();
    const bool initialDecode =
        this->lookupOrDefault<bool>("initialDecode", false);

    forAll(alphas_, phasei)
    {
        // Ensure boundaries are updated
        alphas_[phasei].correctBoundaryConditions();
        rhos_[phasei].correctBoundaryConditions();

        word phaseName = alphas_[phasei].group();
        sharpen_[phasei] = this->subDict(phaseName).lookupOrDefault("sharpen", false);
        fluxScheme_->phases().insert(phaseName);
        alphaRhos_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("alphaRho", phaseName),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                alphas_[phasei]*rhos_[phasei],
                rhos_[phasei].boundaryField().types()
            )
        );
        if (initialDecode && !alphaRhos_[phasei].headerOk())
        {
            alphaRhos_[phasei] = alphas_[phasei]*rhos_[phasei];
            alphaRhos_[phasei].correctBoundaryConditions();
        }

        alphaPhis_.set
        (
            phasei,
            new surfaceScalarField
            (
                IOobject
                (
                    IOobject::groupName("alphaPhi", phaseName),
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0.0)
            )
        );
        alphaRhoPhis_.set
        (
            phasei,
            new surfaceScalarField
            (
                IOobject
                (
                    IOobject::groupName("alphaRhoPhi", phaseName),
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar("0", dimensionSet(1, 0, -1, 0, 0), 0.0)
            )
        );
    }

    if (initialize)
    {
        thermoPtr_->initializeModels();
        this->setModels();

        if (initialDecode && allRead)
        {
            Info<< "Decoding conservative fields"<<endl;
            decode();
        }
        encode();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiphaseCompressibleSystem::~multiphaseCompressibleSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multiphaseCompressibleSystem::update()
{
    decode();

    surfaceScalarField rhoOwn
    (
        surfaceScalarField::New
        (
            reconstruction::ownName("rho"),
            mesh(),
            dimensionedScalar(dimDensity, 0.0)
        )
    );
    surfaceScalarField rhoNei
    (
        surfaceScalarField::New
        (
            reconstruction::neiName("rho"),
            mesh(),
            dimensionedScalar(dimDensity, 0.0)
        )
    );
    PtrList<surfaceScalarField> alphaRhosOwn(alphas_.size());
    PtrList<surfaceScalarField> alphaRhosNei(alphas_.size());
    forAll(alphaRhoPhis_, phasei)
    {
        autoPtr<ReconstructionScheme<scalar>> rhoLimiter
        (
            ReconstructionScheme<scalar>::New
            (
                alphaRhos_[phasei],
                "rho",
                alphaRhos_[phasei].group(),
                true
            )
        );

        alphaRhosOwn.set(phasei, rhoLimiter->interpolateOwn());
        rhoOwn += alphaRhosOwn[phasei];

        alphaRhosNei.set(phasei, rhoLimiter->interpolateNei());
        rhoNei += alphaRhosNei[phasei];
    }

    fluxScheme_->update
    (
        rhoOwn,
        rhoNei,
        U_,
        e_,
        p_,
        speedOfSound()(),
        phi_,
        rhoPhi_,
        rhoUPhi_,
        rhoEPhi_
    );

    tmp<volScalarField> tdivPhi(fvc::div(phi_));
    const volScalarField& divPhi = tdivPhi();

    PtrList<volScalarField> vDots(alphas_.size());

    forAll(alphaRhoPhis_, phasei)
    {
        autoPtr<ReconstructionScheme<scalar>> alphaLimiter
        (
            ReconstructionScheme<scalar>::New
            (
                alphas_[phasei],
                "alpha",
                alphas_[phasei].group(),
                true
            )
        );
        surfaceScalarField alphaOwn(alphaLimiter->interpolateOwn());
        surfaceScalarField alphaNei(alphaLimiter->interpolateNei());

        vDots.set(phasei, alphas_[phasei]*divPhi);
        alphaPhis_[phasei] = fluxScheme_->flux(alphaOwn, alphaNei, phi_);
        alphaRhoPhis_[phasei] = fluxScheme_->flux
        (
            alphaRhosOwn[phasei],
            alphaRhosNei[phasei],
            phi_
        );

        if (mesh().cacheTemporaryObject(alphaOwn.name()))
        {
            mesh().cacheTemporaryObject(alphaOwn);
            mesh().cacheTemporaryObject(alphaNei);
        }
        if (mesh().cacheTemporaryObject(alphaRhosOwn[phasei].name()))
        {
            mesh().cacheTemporaryObject(alphaRhosOwn[phasei]);
            mesh().cacheTemporaryObject(alphaRhosNei[phasei]);
        }
    }

    // Limit alpha flux
    surfaceScalarField phi(phi_);
    this->storeAndBlendDelta(phi);

    UPtrList<const volScalarField> alphas(alphas_.size());
    forAll(alphas_, phasei)
    {
        const volScalarField& alpha = alphas_[phasei];

        volScalarField Sp
        (
            IOobject
            (
                IOobject::groupName("Sp", alpha.group()),
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar(divPhi.dimensions(), 0)
        );

        volScalarField Su
        (
            IOobject::groupName("Su", alpha.group()),
            vDots[phasei]
        );

        {
            const volScalarField& vDot = vDots[phasei];

            forAll(vDot, celli)
            {
                if (vDot[celli] < 0.0 && alpha[celli] > 0.0)
                {
                    Sp[celli] += vDot[celli]*alpha[celli];
                    Su[celli] -= vDot[celli]*alpha[celli];
                }
                else if (vDot[celli] > 0.0 && alpha[celli] < 1.0)
                {
                    Sp[celli] -= vDot[celli]*(1.0 - alpha[celli]);
                }
            }
        }


        forAll(alphas_, phasej)
        {
            if (phasei == phasej) continue;

            const volScalarField& vDot2 = vDots[phasej];
            const volScalarField& alpha2 = alphas_[phasej];

            forAll(vDot2, celli)
            {
                if (vDot2[celli] > 0.0 && alpha2[celli] < 1.0)
                {
                    Sp[celli] -= vDot2[celli]*(1.0 - alpha2[celli]);
                    Su[celli] += vDot2[celli]*alpha[celli];
                }
                else if (vDot2[celli] < 0.0 && alpha2[celli] > 0.0)
                {
                    Sp[celli] += vDot2[celli]*alpha2[celli];
                }
            }
        }
        this->storeAndBlendDelta(Sp);
        this->storeAndBlendDelta(Su);

        surfaceScalarField& alphaPhi = alphaPhis_[phasei];
        this->storeAndBlendDelta(alphaPhi);

        alphas.set(phasei, &alphas_[phasei]);
        volScalarField alphaOld(alphas_[phasei]);
        this->blendOld(alphaOld);
        alphaOld.storeOldTime();

        MULES::limit
        (
            1.0/mesh().time().deltaT().value(),
            geometricOneField(),
            alphaOld,
            phi,
            alphaPhi,
            Sp,
            Su,
            oneField(),
            zeroField(),
            false
        );
        alphaPhi = this->calcAndStoreDelta(alphaPhi);
    }
    MULES::limitSum(alphas, alphaPhis_, phi_);

    thermo_.update();
}


void Foam::multiphaseCompressibleSystem::solve()
{
    // Solve momentum and energy
    compressibleBlastSystem::solve();

    dimensionedScalar dT = rho_.time().deltaT();

    // Divergence of volumetric flux
    volScalarField divPhi(fvc::div(phi_));

    rho_ = Zero;
    forAll(alphaRhos_, phasei)
    {
        volScalarField deltaAlpha
        (
            fvc::div(alphaPhis_[phasei]) - alphas_[phasei]*divPhi
        );
        this->fvTimeInt_->addDeltaSource(alphas_[phasei].name(), deltaAlpha);

        volScalarField deltaAlphaRho(fvc::div(alphaRhoPhis_[phasei]));
        this->fvTimeInt_->addDeltaSource
        (
            alphaRhos_[phasei].name(),
            deltaAlphaRho
        );

        // Blend old values
        this->storeAndBlendOld(alphas_[phasei], false);
        this->storeAndBlendOld(alphaRhos_[phasei]);
        rho_ += alphaRhos_[phasei];

        // Blend deltas
        this->storeAndBlendDelta(deltaAlpha);
        this->storeAndBlendDelta(deltaAlphaRho);

        // Solve volume fraction
        alphas_[phasei] -= dT*deltaAlpha;
        alphas_[phasei].maxMin(0.0, 1.0);
        alphas_[phasei].correctBoundaryConditions();

        // Solve phase mass transport
        alphaRhos_[phasei].storePrevIter();
        alphaRhos_[phasei] -= dT*deltaAlphaRho;
        alphaRhos_[phasei].correctBoundaryConditions();

        if (transportPhaseDensity_)
        {
            volScalarField deltaRho
            (
                IOobject::groupName("deltaRho", rhos_[phasei].group()),
                fvc::div(fluxScheme_->flux(rhos_[phasei], phi_))
              - rhos_[phasei]*divPhi
            );

            this->storeAndBlendDelta(deltaRho);
            this->storeAndBlendOld(rhos_[phasei], false);

            //- Solve volume fraction
            rhos_[phasei] -= dT*deltaRho;
            rhos_[phasei].correctBoundaryConditions();
        }
    }

    // Store "old" total density
    rho_.storePrevIter();

    //- Compute new density
    rho_ = Zero;
    forAll(alphas_, phasei)
    {
        rho_ += alphaRhos_[phasei];
    }

    // Solve thermo
    thermoPtr_->solve();
}


void Foam::multiphaseCompressibleSystem::postUpdate()
{
    this->decode();

    // Solve volume fraction and phase mass transports
    bool needUpdate = false;
    forAll(alphas_, phasei)
    {
        if (needSolve(alphas_[phasei].name()))
        {
            needUpdate = true;

            fvScalarMatrix alphaEqn
            (
                fvm::ddt(alphas_[phasei]) - fvc::ddt(alphas_[phasei])
             ==
                models().source(alphas_[phasei])
            );
            constraints().constrain(alphaEqn);
            alphaEqn.solve();
            constraints().constrain(alphas_[phasei]);
        }
    }
    if (needUpdate)
    {
        calcAlphas();
    }

    // Solve phase 1 mass
    rho_.storePrevIter();
    rho_ = dimensionedScalar(dimDensity, 0.0);
    forAll(rhos_, phasei)
    {
        volScalarField& rho(rhos_[phasei]);
        if (needSolve(rho.name()))
        {
            const volScalarField& alpha(alphas_[phasei]);
            dimensionedScalar rAlpha(thermo_.thermo(phasei).residualAlpha());
            fvScalarMatrix alphaRhoEqn
            (
                fvm::ddt(alpha, rho) - fvc::ddt(alphaRhos_[phasei])
              + fvm::ddt(rAlpha, rho) - fvc::ddt(rAlpha, rho)
            ==
                models().source(alpha, rho)
            );
            constraints().constrain(alphaRhoEqn);
            alphaRhoEqn.solve();
            constraints().constrain(rho);

            alphaRhos_[phasei] = alpha*rho;
        }
        rho_ += alphaRhos_[phasei];
    }

    compressibleBlastSystem::postUpdate();
}


void Foam::multiphaseCompressibleSystem::calcAlphas()
{
    // find largest volume fraction and set to 1-sum
    SortableList<scalar> alphas(alphas_.size());
    forAll(rho_, celli)
    {
        scalar sumAlpha = 0.0;
        forAll(alphas_, phasei)
        {
            alphas[phasei] = alphas_[phasei][celli];
            sumAlpha += alphas_[phasei][celli];
        }

        alphas.reverseSort();

        const labelList& indices = alphas.indices();
        const label fixedPhase = indices[0];
        sumAlpha -= alphas[0];
        if (sumAlpha > 1)
        {
            for (label i = 1; i < alphas.size(); i++)
            {
                const label phasei = indices[i];
                alphas_[phasei][celli] /= sumAlpha;
            }
            alphas_[fixedPhase][celli] = 0.0;
        }
        else
        {
            alphas_[fixedPhase][celli] = 1.0 - sumAlpha;
        }
    }
    forAll(alphas_, phasei)
    {
        alphas_[phasei].correctBoundaryConditions();
    }
}


void Foam::multiphaseCompressibleSystem::decode()
{
    calcAlphas();

    // Calculate densities
    rho_ = Zero;

    forAll(alphas_, phasei)
    {
        volScalarField& alpha = alphas_[phasei];
        volScalarField& rho = rhos_[phasei];
        volScalarField& alphaRho = alphaRhos_[phasei];
        const scalar rAlpha = thermo_.thermo(phasei).residualAlpha().value();

        alphaRho.max(0);
        if (transportPhaseDensity_)
        {
            // Only update cells that have a valid volume fraction
            // other cell densities are handled by transport of density
            forAll(alpha, celli)
            {
                const scalar alphai = alpha[celli];
                if (alphai > rAlpha)
                {
                    rho[celli] = alphaRho[celli]/alphai;
                }
            }
        }
        else
        {
            rho.ref() = alphaRho()/max(alpha(), rAlpha);
        }

        rho.correctBoundaryConditions();
        alphaRho.boundaryFieldRef() =
            alpha.boundaryField()*rho.boundaryField();

        rho_ += alphaRhos_[phasei];
    }

    compressibleBlastSystem::decode();
}


void Foam::multiphaseCompressibleSystem::encode()
{
    rho_ = Zero;
    forAll(alphas_, phasei)
    {
        alphaRhos_[phasei] = alphas_[phasei]*rhos_[phasei];
        rho_ += alphaRhos_[phasei];
    }
    compressibleBlastSystem::encode();
}

// ************************************************************************* //
