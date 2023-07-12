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
    const fvMesh& mesh
)
:
    compressibleBlastSystem(mesh, multiphaseFluidBlastThermo::typeName),
    thermo_(dynamicCast<multiphaseFluidBlastThermo>(thermoPtr_())),
    alphas_(thermo_.volumeFractions()),
    rhos_(thermo_.rhos()),
    alphaRhos_(alphas_.size()),
    alphaPhis_(alphas_.size()),
    alphaRhoPhis_(alphas_.size()),
    transportPhaseDensity_(this->lookupOrDefault("transportPhaseDensity", false))
{
    this->fluxScheme_ = fluxScheme::NewMulti(phi_);

    forAll(alphas_, phasei)
    {
        // Ensure boundaries are updated
        alphas_[phasei].correctBoundaryConditions();
        alphas_[phasei].oldTime();

        rhos_[phasei].correctBoundaryConditions();
        rhos_[phasei].oldTime();

        word phaseName = alphas_[phasei].group();
        alphaRhos_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("alphaRho", phaseName),
                    mesh.time().timeName(),
                    mesh
                ),
                alphas_[phasei]*rhos_[phasei],
                rhos_[phasei].boundaryField().types()
            )
        );
        alphaRhos_[phasei].oldTime();

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

    thermoPtr_->initializeModels();
    this->setModels();

    encode();
}


Foam::multiphaseCompressibleSystem::multiphaseCompressibleSystem
(
    const fvMesh& mesh,
    const bool
)
:
    compressibleBlastSystem(mesh, multiphaseFluidBlastThermo::typeName),
    thermo_(dynamicCast<multiphaseFluidBlastThermo>(thermoPtr_())),
    alphas_(thermo_.volumeFractions()),
    rhos_(thermo_.rhos()),
    alphaRhos_(alphas_.size()),
    alphaPhis_(alphas_.size()),
    alphaRhoPhis_(alphas_.size()),
    transportPhaseDensity_(this->lookupOrDefault("transportPhaseDensity", false))
{
    this->fluxScheme_ = fluxScheme::NewMulti(phi_);

    forAll(alphas_, phasei)
    {
        // Ensure boundaries are updated
        alphas_[phasei].correctBoundaryConditions();
        alphas_[phasei].oldTime();

        rhos_[phasei].correctBoundaryConditions();
        rhos_[phasei].oldTime();

        word phaseName = alphas_[phasei].group();
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
                    mesh
                ),
                alphas_[phasei]*rhos_[phasei],
                rhos_[phasei].boundaryField().types()
            )
        );
        alphaRhos_[phasei].oldTime();

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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiphaseCompressibleSystem::~multiphaseCompressibleSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multiphaseCompressibleSystem::update()
{
    decode();
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

        alphaPhis_[phasei] = fluxScheme_->flux(talphaOwn(), talphaNei(), phi_);

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

        tmp<surfaceScalarField> talphaRhoOwn
        (
            surfaceScalarField::New
            (
                rhoLimiter->ownName(alphaRhos_[phasei].name()),
                talphaOwn*trhoOwn()
            )
        );
        tmp<surfaceScalarField> talphaRhoNei
        (
            surfaceScalarField::New
            (
                rhoLimiter->neiName(alphaRhos_[phasei].name()),
                talphaNei*trhoNei()
            )
        );

        alphaRhoPhis_[phasei] =
            fluxScheme_->flux(talphaRhoOwn(), talphaRhoNei(), phi_);
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
        this->storeAndBlendDelta(deltaAlpha);
        this->storeAndBlendOld(alphas_[phasei], false);

        volScalarField deltaAlphaRho(fvc::div(alphaRhoPhis_[phasei]));
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
    // find largest volume fraction and set to 1-sum
    forAll(rho_, celli)
    {
        SortableList<scalar> alphas(alphas_.size());
        forAll(alphas_, phasei)
        {
            alphas_[phasei][celli] = max(min(alphas_[phasei][celli], 1.0), 0.0);
            alphas[phasei] = alphas_[phasei][celli];
        }
        alphas.reverseSort();

        const label fixedPhase = alphas.indices()[0];

        scalar sumAlpha = 0.0;
        for (label phasei = 1; phasei < alphas.size(); phasei++)
        {
            sumAlpha += alphas[phasei];
        }
        alphas_[fixedPhase][celli] = 1.0 - sumAlpha;
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

        alphaRho.boundaryFieldRef() = alpha.boundaryField()*rho.boundaryField();

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
