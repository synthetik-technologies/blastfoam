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

#include "interfaceTwoPhaseCompressibleSystem.H"
#include "meshSizeObject.H"
#include "fvcCellReduce.H"
#include "upwind.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(interfaceTwoPhaseCompressibleSystem, 0);
    addToRunTimeSelectionTable
    (
        compressibleSystem,
        interfaceTwoPhaseCompressibleSystem,
        twoPhase
    );
}


// * * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceTwoPhaseCompressibleSystem::interfaceTwoPhaseCompressibleSystem
(
    const fvMesh& mesh
)
:
    compressibleBlastSystem(2, mesh),
    thermo_
    (
        refCast<twoPhaseFluidBlastThermo>(thermoPtr_())
    ),
    alpha1_(thermo_.alpha1()),
    rho1_(thermo_.thermo(0).rho()),
    alpha2_(thermo_.alpha2()),
    rho2_(thermo_.thermo(1).rho()),
    levelSet_(alpha1_, this->subDict(alpha1_.group())),
    redistance_(this->lookupOrDefault("redistance", false)),
    redistanceInterval_
    (
        redistance_
      ? this->lookup<label>("redistanceInterval")
      : labelMax
    ),
    alphaRho1_
    (
        IOobject
        (
            IOobject::groupName("alphaRho", rho1_.group()),
            mesh.time().timeName(),
            mesh
        ),
        alpha1_*rho1_
    ),
    alphaRho2_
    (
        IOobject
        (
            IOobject::groupName("alphaRho", rho2_.group()),
            mesh.time().timeName(),
            mesh
        ),
        alpha2_*rho2_
    ),
    alphaPhi1_
    (
        IOobject
        (
            IOobject::groupName("alphaPhi", alphaRho1_.group()),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0.0)
    ),
    alphaPhi2_
    (
        IOobject
        (
            IOobject::groupName("alphaPhi", alphaRho2_.group()),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0.0)
    ),
    alphaRhoPhi1_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoPhi", alphaRho1_.group()),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(1, 0, -1, 0, 0), 0.0)
    ),
    alphaRhoPhi2_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoPhi", alphaRho2_.group()),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(1, 0, -1, 0, 0), 0.0)
    ),
    sigmaPtr_(surfaceTensionModel::New(*this, mesh))
{
    alpha1_ = levelSet_.alpha();
    alpha2_ = 1.0 - alpha1_;
    thermo_.initializeFields();
    thermo_.correct();
    this->fluxScheme_ = fluxScheme::NewMulti(mesh);

    rho_ = alphaRho1_ + alphaRho2_;

    thermo_.initializeModels();
    this->setModels();

    encode();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceTwoPhaseCompressibleSystem::~interfaceTwoPhaseCompressibleSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::interfaceTwoPhaseCompressibleSystem::update()
{
    // Correct the levelSet variables
    if
    (
        this->step() == 0
     && (mesh().time().timeIndex() % redistanceInterval_ == 0)
    )
    {
        levelSet_.redistance();
    }
    levelSet_.correct();

    decode();
    thermo_.update();

    fluxScheme_->update
    (
        alpha1_,
        rho1_,
        rho2_,
        U_,
        e_,
        p_,
        speedOfSound()(),
        phi_,
        alphaPhi1_,
        alphaRhoPhi1_,
        alphaRhoPhi2_,
        rhoPhi_,
        rhoUPhi_,
        rhoEPhi_
    );

    alphaPhi2_ = phi_ - alphaPhi1_;
}


void Foam::interfaceTwoPhaseCompressibleSystem::solve()
{
    dimensionedScalar dT = rho_.time().deltaT();

    volScalarField& levelSet(levelSet_.levelSet());

    //- Update changes in volume fraction and phase mass
    volScalarField deltaLevelSet
    (
        fvc::div(fluxScheme_->interpolate(levelSet)*phi_)
      - levelSet*fvc::div(phi_)
    );
    volScalarField deltaAlphaRho1(fvc::div(alphaRhoPhi1_));
    volScalarField deltaAlphaRho2(fvc::div(alphaRhoPhi2_));

    this->storeAndBlendDelta(deltaLevelSet);
    this->storeAndBlendDelta(deltaAlphaRho1);
    this->storeAndBlendDelta(deltaAlphaRho2);

    // Volume fraction and levelSet fields are not scaled by change in volume
    // because it is not a conservative quantity
    this->storeAndBlendOld(levelSet, false);

    this->storeAndBlendOld(alphaRho1_);
    alphaRho1_.storePrevIter();

    this->storeAndBlendOld(alphaRho2_);
    alphaRho2_.storePrevIter();

    rho_ = alphaRho1_ + alphaRho2_;
    rho_.storePrevIter();

    levelSet -= dT*deltaLevelSet;
    levelSet.correctBoundaryConditions();

    alphaRho1_ -= dT*deltaAlphaRho1;
    alphaRho1_.correctBoundaryConditions();

    alphaRho2_ -= dT*deltaAlphaRho2;
    alphaRho2_.correctBoundaryConditions();

    rho_ = alphaRho1_ + alphaRho2_;

    thermo_.solve();

    //- Calculate deltas for momentum and energy
    volVectorField sTF
    (
        fvc::reconstruct
        (
            fvc::interpolate(sigmaPtr_->sigma()*levelSet_.K())
           *fvc::snGrad(alpha1_)
           *rho_.mesh().magSf()
        )
    );
    volVectorField deltaRhoU
    (
        "deltaRhoU",
        fvc::div(rhoUPhi_) - g_*rho_ + sTF
    );

    volScalarField deltaRhoE
    (
        "deltaRhoE",
        fvc::div(rhoEPhi_)
      - ESource()
      - (rhoU_ & g_)
      + (U_ & sTF)
    );

    //- Store old values
    this->storeAndBlendOld(rhoU_);
    this->storeAndBlendOld(rhoE_);

    //- Store changed in momentum and energy
    this->storeAndBlendDelta(deltaRhoU);
    this->storeAndBlendDelta(deltaRhoE);

    //- Solve for momentum and energy
    rhoU_ -= cmptMultiply(dT*deltaRhoU, solutionDs_);
    rhoE_ -= dT*deltaRhoE;
}


void Foam::interfaceTwoPhaseCompressibleSystem::postUpdate()
{
    this->decode();

    // Solve volume fraction
    if (needSolve(alpha1_.name()))
    {
        fvScalarMatrix alphaEqn
        (
            fvm::ddt(alpha1_) - fvc::ddt(alpha1_)
        ==
            models().source(alpha1_)
        );
        constraints().constrain(alphaEqn);
        alphaEqn.solve();
        constraints().constrain(alpha1_);

        alpha1_.maxMin(0.0, 1.0);
        alpha2_ = 1.0 - alpha1_;
    }
    // Solve phase 1 mass
    if (needSolve(rho1_.name()))
    {
        fvScalarMatrix alphaRho1Eqn
        (
            fvm::ddt(alpha1_, rho1_) - fvc::ddt(alphaRho1_)
          + fvm::ddt(thermoPtr_->residualAlpha(), rho1_)
          - fvc::ddt(thermoPtr_->residualAlpha(), rho1_)
        ==
            models().source(alpha1_, rho1_)
        );
        constraints().constrain(alphaRho1Eqn);
        alphaRho1Eqn.solve();
        constraints().constrain(rho1_);
    }
    // Solve phase 2 mass
    if (needSolve(rho2_.name()))
    {
        fvScalarMatrix alphaRho2Eqn
        (
            fvm::ddt(alpha2_, rho2_) - fvc::ddt(alphaRho2_)
          + fvm::ddt(thermoPtr_->residualAlpha(), rho2_)
          - fvc::ddt(thermoPtr_->residualAlpha(), rho2_)
        ==
            models().source(alpha2_, rho2_)
        );
        constraints().constrain(alphaRho2Eqn);
        alphaRho2Eqn.solve();
        constraints().constrain(rho2_);
    }

    // Update phase masses
    alphaRho1_ = alpha1_*rho1_;
    alphaRho2_ = alpha2_*rho2_;

    rho_.storePrevIter();
    rho_ = alphaRho1_ + alphaRho2_;

    compressibleBlastSystem::postUpdate();
}


void Foam::interfaceTwoPhaseCompressibleSystem::decode()
{
    alpha1_ = levelSet_.alpha();
    alpha1_.correctBoundaryConditions();
    alpha2_ = 1.0 - alpha1_;
    alpha1_.correctBoundaryConditions();


    // Stabilized fields
    volScalarField alpha1
    (
        max(alpha1_, thermo_.thermo(0).residualAlpha())
    );
    volScalarField alpha2
    (
        max(alpha2_, thermo_.thermo(1).residualAlpha())
    );

    volScalarField a1Pos(alpha1_ - thermo_.thermo(0).residualAlpha());
    volScalarField a2Pos(1.0 - a1Pos);
    alphaRho1_.max(0);
    rho1_.ref() = (alphaRho1_()/alpha1())*a1Pos() + a2Pos()*rho1_();
    rho1_.max(thermo_.thermo(0).residualRho());
    rho1_.correctBoundaryConditions();
    alphaRho1_.boundaryFieldRef() = alpha1_.boundaryField()*rho1_.boundaryField();

    alphaRho2_.max(0);
    rho2_.ref() = (alphaRho2_()/alpha2())*a2Pos() + a1Pos()*rho2_();
    rho2_.max(thermo_.thermo(1).residualRho());
    rho2_.correctBoundaryConditions();
    alphaRho2_.boundaryFieldRef() = alpha2_.boundaryField()*rho2_.boundaryField();

    rho_ == alphaRho1_ + alphaRho2_;

    compressibleBlastSystem::decode();
}


void Foam::interfaceTwoPhaseCompressibleSystem::encode()
{
    alphaRho1_ = alpha1_*rho1_;
    alphaRho2_ = alpha2_*rho2_;
    rho_ = alphaRho1_ + alphaRho2_;
    compressibleBlastSystem::encode();
}

// ************************************************************************* //
