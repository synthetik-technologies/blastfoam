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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::singlePhaseCompressibleSystem::singlePhaseCompressibleSystem
(
    const fvMesh& mesh
)
:
    compressibleBlastSystem(1, mesh)
{
    this->fluxScheme_ = fluxScheme::NewSingle(mesh);

    thermoPtr_->initializeModels();
    this->setModels();
    encode();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::singlePhaseCompressibleSystem::~singlePhaseCompressibleSystem()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::singlePhaseCompressibleSystem::solve()
{
    volScalarField deltaRho("deltaRho", fvc::div(rhoPhi_));
    this->storeAndBlendDelta(deltaRho);

    dimensionedScalar dT = rho_.time().deltaT();
    this->storeAndBlendOld(rho_);
    rho_.storePrevIter();

    rho_ -= dT*deltaRho;
    rho_.correctBoundaryConditions();

    thermoPtr_->solve();

    compressibleBlastSystem::solve();
}


void Foam::singlePhaseCompressibleSystem::postUpdate()
{
    this->decode();

    // Solve mass
    if (needSolve(rho_.name()))
    {
        fvScalarMatrix rhoEqn
        (
            fvm::ddt(rho_) - fvc::ddt(rho_)
         ==
            models().source(rho_)
        );
        constraints().constrain(rhoEqn);
        rhoEqn.solve();
        constraints().constrain(rho_);
    }

    compressibleBlastSystem::postUpdate();
}

// ************************************************************************* //
