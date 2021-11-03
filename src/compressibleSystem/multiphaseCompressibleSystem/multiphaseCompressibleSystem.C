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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiphaseCompressibleSystem, 0);
    addToRunTimeSelectionTable
    (
        compressibleSystem,
        multiphaseCompressibleSystem,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseCompressibleSystem::multiphaseCompressibleSystem
(
    const fvMesh& mesh
)
:
    compressibleBlastSystem(3, mesh),
    thermo_(dynamicCast<multiphaseFluidBlastThermo>(thermoPtr_())),
    alphas_(thermo_.volumeFractions()),
    rhos_(thermo_.rhos()),
    alphaRhos_(alphas_.size()),
    alphaPhis_(alphas_.size()),
    alphaRhoPhis_(alphas_.size())
{
    forAll(alphas_, phasei)
    {
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


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiphaseCompressibleSystem::~multiphaseCompressibleSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multiphaseCompressibleSystem::update()
{
    decode();
    fluxScheme_->update
    (
        alphas_,
        rhos_,
        U_,
        e_,
        p_,
        speedOfSound()(),
        phi_,
        alphaPhis_,
        alphaRhoPhis_,
        rhoPhi_,
        rhoUPhi_,
        rhoEPhi_
    );
    thermo_.update();
}


void Foam::multiphaseCompressibleSystem::solve()
{
    dimensionedScalar dT = rho_.time().deltaT();
    rho_ = dimensionedScalar("0", dimDensity, 0.0);
    forAll(alphas_, phasei)
    {
        volScalarField deltaAlpha
        (
            fvc::div(alphaPhis_[phasei]) - alphas_[phasei]*fvc::div(phi_)
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
                fvm::ddt(alpha, rho) - fvc::ddt(alpha, rho)
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
    // Calculate densities
    rho_ = dimensionedScalar("0", dimDensity, 0.0);

    forAll(alphas_, phasei)
    {
        alphaRhos_[phasei].max(0);
        rhos_[phasei] =
            alphaRhos_[phasei]
           /max(alphas_[phasei], thermo_.thermo(phasei).residualAlpha());

        rhos_[phasei].correctBoundaryConditions();

        alphaRhos_[phasei] = alphas_[phasei]*rhos_[phasei];

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
