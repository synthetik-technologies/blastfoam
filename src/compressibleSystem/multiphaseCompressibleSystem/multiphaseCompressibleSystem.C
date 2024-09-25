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

    PtrList<surfaceScalarField> alphasOwn(alphas_.size());
    PtrList<surfaceScalarField> alphasNei(alphas_.size());
    PtrList<surfaceScalarField> alphaRhosOwn(alphas_.size());
    PtrList<surfaceScalarField> alphaRhosNei(alphas_.size());
    tmp<surfaceScalarField> trhoOwn
    (
        surfaceScalarField::New
        (
            reconstruction::ownName(rho_.name()),
            mesh(),
            dimensionedScalar(dimDensity, 0.0)
        )
    );
    surfaceScalarField& rhoOwn = trhoOwn.ref();
    tmp<surfaceScalarField> trhoNei
    (
        surfaceScalarField::New
        (
            reconstruction::neiName(rho_.name()),
            mesh(),
            dimensionedScalar(dimDensity, 0.0)
        )
    );
    surfaceScalarField& rhoNei = trhoNei.ref();

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
        tmp<surfaceScalarField> talphaOwn, talphaNei;
        alphaLimiter->interpolateOwnNei(talphaOwn, talphaNei);

        autoPtr<ReconstructionScheme<scalar>> rhoLimiter
        (
            ReconstructionScheme<scalar>::New
            (
                rhos_[phasei],
                "rho",
                rhos_[phasei].group(),
                true
            )
        );
        tmp<surfaceScalarField> trhoOwn, trhoNei;
        rhoLimiter->interpolateOwnNei(trhoOwn, trhoNei);
        fluxScheme::correctPhaseFields
        (
            alphas_[phasei],
            trhoOwn.ref(), trhoNei.ref(),
            thermo_.thermo(phasei).residualAlpha().value()
        );

        alphaRhosOwn.set
        (
            phasei,
            surfaceScalarField::New
            (
                rhoLimiter->ownName(alphaRhos_[phasei].name()),
                talphaOwn()*trhoOwn()
            )
        );
        alphaRhosNei.set
        (
            phasei,
            surfaceScalarField::New
            (
                rhoLimiter->neiName(alphaRhos_[phasei].name()),
                talphaNei()*trhoNei()
            )
        );
        alphasOwn.set(phasei, talphaOwn.ptr());
        alphasNei.set(phasei, talphaNei.ptr());
        rhoOwn += alphaRhosOwn[phasei];
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

    forAll(alphaRhoPhis_, phasei)
    {
        alphaPhis_[phasei] =
            fluxScheme_->flux(alphasOwn[phasei], alphasNei[phasei], phi_);
        alphaRhoPhis_[phasei] =
            fluxScheme_->flux(alphaRhosOwn[phasei], alphaRhosNei[phasei], phi_);
    }
    thermo_.update();
}


void Foam::multiphaseCompressibleSystem::solve()
{
    dimensionedScalar dT = rho_.time().deltaT();
    rho_ = dimensionedScalar("0", dimDensity, 0.0);

    volScalarField divPhi(fvc::div(phi_));
    forAll(alphas_, phasei)
    {
        volScalarField deltaAlpha
        (
            fvc::div(alphaPhis_[phasei]) - alphas_[phasei]*divPhi
        );
        this->fvTimeInt_->addDeltaSource(alphas_[phasei].name(), deltaAlpha);
        this->storeAndBlendDelta(deltaAlpha);
        this->storeAndBlendOld(alphas_[phasei], false);

        volScalarField deltaAlphaRho(fvc::div(alphaRhoPhis_[phasei]));
        this->fvTimeInt_->addDeltaSource
        (
            alphaRhos_[phasei].name(),
            deltaAlphaRho
        );
        this->storeAndBlendDelta(deltaAlphaRho);
        this->storeAndBlendOld(alphaRhos_[phasei]);
        rho_ += alphaRhos_[phasei];

        //- Solve volume fraction
        alphas_[phasei] -= dT*deltaAlpha;
        alphas_[phasei].correctBoundaryConditions();


        //- Solve phasic mass transport
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

    //- Store "old" total density
    rho_.storePrevIter();

    //- Compute new density
    rho_ = dimensionedScalar("0", dimDensity, 0.0);
    forAll(alphas_, phasei)
    {
        rho_ += alphaRhos_[phasei];
    }

    thermoPtr_->solve();

    compressibleBlastSystem::solve();
}


void Foam::multiphaseCompressibleSystem::postUpdate()
{
    // Solve volume fraction and phase mass transports
    bool needUpdate = false;
    surfaceScalarField nMagUf(fvc::interpolate(mag(U_)));
    nMagUf.dimensions().reset(sqr(dimLength)/dimTime);
    forAll(alphas_, phasei)
    {
        if (needSolve(alphas_[phasei].name()) || sharpen_[phasei])
        {
            needUpdate = true;

            if (sharpen_[phasei])
            {
                const dictionary& controls = mesh().solverDict(alphas_[phasei].name());
                const label nCorrectors
                (
                    controls.lookupOrDefault<label>("nCorrectors", 10)
                );
                const scalar sharpenTol(controls.lookup<scalar>("sharpenTol"));
                label iter = 0;
                scalar initialResidual = 0;
                do
                {
                    fvScalarMatrix alphaEqn
                    (
                        fvm::ddt(alphas_[phasei])
                      - fvc::ddt(alphas_[phasei])
                      + fvm::laplacian(nMagUf, alphas_[phasei])
                    );
                    initialResidual = alphaEqn.solve().initialResidual();
                } while (++iter < nCorrectors && initialResidual > sharpenTol);
            }

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

    // Update primitive variables
    this->decode();

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
    // // find largest volume fraction and set to 1-sum
    // forAll(rho_, celli)
    // {
    //     // SortableList<scalar> alphas(alphas_.size());
    //     scalar sumAlpha = 0.0;
    //     forAll(alphas_, phasei)
    //     {
    //         alphas_[phasei][celli] = max(alphas_[phasei][celli], 0.0);
    //         sumAlpha += sqr(alphas_[phasei][celli]);
    //         // alphas_[phasei][celli] = max(min(alphas_[phasei][celli], 1.0), 0.0);
    //         // alphas[phasei] = alphas_[phasei][celli];
    //     }
    //     sumAlpha = sqrt(sumAlpha);
    //     forAll(alphas_, phasei)
    //     {
    //         alphas_[phasei][celli] /= sumAlpha;
    //     }
    //
    //     // alphas.reverseSort();
    //
    //     // const label fixedPhase = alphas.indices()[0];
    //
    //     // for (label phasei = 1; phasei < alphas.size(); phasei++)
    //     // {
    //     //     sumAlpha += alphas[phasei];
    //     // }
    //     // alphas_[fixedPhase][celli] = 1.0 - sumAlpha;
    // }
    forAll(rho_, celli)
    {
        scalar sumAlpha = 1.0;
        forAll(alphas_, phasei)
        {
            alphas_[phasei][celli] = min(max(alphas_[phasei][celli], 0.0), sumAlpha);
            sumAlpha = max(sumAlpha - alphas_[phasei][celli], 0);
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
    rho_ = dimensionedScalar("0", dimDensity, 0.0);

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
        alphaRho.correctBoundaryConditions();
        alphaRho.boundaryFieldRef() =
            alpha.boundaryField()*rho.boundaryField();

        rho_ += alphaRhos_[phasei];
    }

    compressibleBlastSystem::decode();
}


void Foam::multiphaseCompressibleSystem::encode()
{
    rho_ = dimensionedScalar("0", dimDensity, 0.0);
    forAll(alphas_, phasei)
    {
        alphaRhos_[phasei] = alphas_[phasei]*rhos_[phasei];
        rho_ += alphaRhos_[phasei];
    }
    compressibleBlastSystem::encode();
}

// ************************************************************************* //
