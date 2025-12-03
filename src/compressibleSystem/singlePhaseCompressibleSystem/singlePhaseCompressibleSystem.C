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

#include "singlePhaseCompressibleSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(singlePhaseCompressibleSystem, 0);
    addToRunTimeSelectionTable
    (
        compressibleSystem,
        singlePhaseCompressibleSystem,
        singlePhase
    );
}

// * * * * * * * * * * * * Protected Members Functions * * * * * * * * * * * //

void Foam::singlePhaseCompressibleSystem::solveMass()
{
    volScalarField& rho = this->rhoEff();
    dimensionedScalar dT = rho.time().deltaT();

    volScalarField deltaRho("deltaRho", fvc::div(rhoPhi_));
    this->fvTimeInt_->addDeltaSource(rho_.name(), deltaRho);

    if (this->LTS())
    {
        deltaRho /= corDeltaT();
    }

    this->storeAndBlendDelta(deltaRho);
    this->storeAndBlendOld(rho);

    rho.storePrevIter();

    rho -= dT*deltaRho;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::singlePhaseCompressibleSystem::singlePhaseCompressibleSystem
(
    const dictionary& dict,
    const fvMesh& mesh,
    const bool initialize
)
:
    compressibleBlastSystem(dict, mesh, word::null)
{
    this->fluxScheme_ = fluxScheme::NewSingle(phi_);

    if (initialize)
    {
        thermoPtr_->initializeModels();
        this->setModels();
        encode();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::singlePhaseCompressibleSystem::~singlePhaseCompressibleSystem()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::singlePhaseCompressibleSystem::decode()
{
    this->rhoEff().correctBoundaryConditions();
    compressibleBlastSystem::decode();
}


void Foam::singlePhaseCompressibleSystem::postImplicit()
{
    this->decode();

    // Solve mass
    volScalarField& rho = this->rhoEff();
    rho.storePrevIter();
    if (needSolve(rho.name()) || rhoSource_.valid())
    {
        fvScalarMatrix rhoEqn
        (
            fvm::ddt(rho) - rhoAdvection_()
         ==
            models().source(rho)
        );

        if (rhoSource_.valid())
        {
            rhoEqn -= rhoSource_;
        }

        constraints().constrain(rhoEqn);
        rhoEqn.solve(rho_.name());
        constraints().constrain(rho);
    }

    compressibleBlastSystem::postImplicit();
}


void Foam::singlePhaseCompressibleSystem::storeExplicit()
{
    compressibleBlastSystem::storeExplicit();

    if (needSolve(rhoEff().name()) || rhoSource_.valid())
    {
        rhoAdvection_ = fvc::ddt(rhoEff());
    }
}


void Foam::singlePhaseCompressibleSystem::clear()
{
    compressibleBlastSystem::clear();

    rhoAdvection_.clear();
}

// ************************************************************************* //
