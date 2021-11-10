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


Foam::tmp<Foam::volScalarField>
Foam::interfaceTwoPhaseCompressibleSystem::H
(
    const volScalarField& levelSet
) const
{
    const volScalarField& h = meshSizeObject::New(mesh()).dx();
    return 0.5*(1.0 + tanh(levelSet_/(h*epsilon_)));


//     scalar pi = Foam::constant::mathematical::pi;
//     volScalarField e(h*epsilon_);
//     tmp<volScalarField> talpha
//     (
//         0.5
//       + levelSet/2.0*e
//       + 1.0/2.0/pi*sin(pi*levelSet/e)
//     );
//     volScalarField& alpha = talpha.ref();
//     forAll(alpha, celli)
//     {
//         if (levelSet[celli] > e[celli])
//         {
//             alpha[celli] = 1.0;
//         }
//         else if (levelSet[celli] < -e[celli]
//         {
//             alpha[celli] = 0.0;
//         }
//     }
//     return talpha;
}


Foam::tmp<Foam::volScalarField>
Foam::interfaceTwoPhaseCompressibleSystem::Hinv
(
    const volScalarField& alpha
) const
{
    const volScalarField& h = meshSizeObject::New(mesh()).dx();
//     volScalarField e(h*epsilon_);
    return h*epsilon_*atanh(2.0*min(0.999, max(0.001, alpha)) - 1.0);
}


Foam::tmp<Foam::volScalarField>
Foam::interfaceTwoPhaseCompressibleSystem::sharpen
(
    const volScalarField& alpha
) const
{
    return
        min
        (
            max
            (
                0.5*((alpha - 0.5)/(0.5 - epsilon_) + 1.0),
                0.0
            ),
            1.0
        );
}


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
    epsilon_("epsilon", dimless, *this),
    levelSet_
    (
        IOobject
        (
            IOobject::groupName("levelSet", rho1_.group()),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        Hinv(alpha1_)
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
    alphaPhi_
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
    interfaceProperties_
    (
        alpha1_,
        alpha2_,
        U_,
        *this
    )
{
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
    decode();

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
        alphaPhi_,
        alphaRhoPhi1_,
        alphaRhoPhi2_,
        rhoPhi_,
        rhoUPhi_,
        rhoEPhi_
    );
    interfaceProperties_.correct();
    thermo_.update();
}


void Foam::interfaceTwoPhaseCompressibleSystem::solve()
{

    //- Update changes in volume fraction and phase mass
    volScalarField deltaAlpha
    (
        fvc::div(alphaPhi_)
      - alpha1_*fvc::div(phi_)
    );
    volScalarField deltaAlphaRho1(fvc::div(alphaRhoPhi1_));
    volScalarField deltaAlphaRho2(fvc::div(alphaRhoPhi2_));

    this->storeAndBlendDelta(deltaAlpha);
    this->storeAndBlendDelta(deltaAlphaRho1);
    this->storeAndBlendDelta(deltaAlphaRho2);


    dimensionedScalar dT = rho_.time().deltaT();

    // Update conservative quantities
    this->storeAndBlendOld(alpha1_);
    this->storeAndBlendOld(alphaRho1_);
    this->storeAndBlendOld(alphaRho2_);
    rho_ = alphaRho1_ + alphaRho2_;
    rho_.storePrevIter();

    alpha1_ -= dT*deltaAlpha;
    alpha1_.correctBoundaryConditions();

    alphaRho1_.storePrevIter();
    alphaRho1_ -= dT*deltaAlphaRho1;
    alphaRho1_.correctBoundaryConditions();

    alphaRho2_.storePrevIter();
    alphaRho2_ -= dT*deltaAlphaRho2;
    alphaRho2_.correctBoundaryConditions();

    rho_ = alphaRho1_ + alphaRho2_;

    // Solve level set transport and correct volume fraction
//     const volScalarField& h = meshSizeObject::New(mesh()).dx();
//     for (label i = 0; i < 2; i++)
//     {
// //         volScalarField Gamma
// //         (
// //             2.0*h*epsilon_/sqr(max(2.0*alpha1_ - 1.0, thermo_.residualAlpha()))
// //            *(
// //                 1.0
// //                /max
// //                 (
// //                     1.0 - sqr(2.0*alpha1_ - 1.0),
// //                     thermo_.residualAlpha()
// //                 ) - 1.0
// //              )
// //         );
// //         Gamma.max(small);
// //         surfaceScalarField levelSetf
// //         (
// //             fluxScheme_->interpolate(levelSet_, "levelSet")
// //         );
// //         volScalarField deltaLevelSet
// //         (
// //             "deltaLevelSet",
// //             fvc::div(levelSetf*phi_) - levelSet_*fvc::div(phi_)
// //         );
// //
// //         this->storeAndBlendDelta(deltaLevelSet);
// //         this->storeAndBlendOld(levelSet_, false);
// //         levelSet_ -= dT*deltaLevelSet;
// //
// //         levelSet_ += sign(levelSet_)*(1.0 - mag(fvc::grad(levelSetf)));
// //         volScalarField magGradAlpha
// //         (
// //             mag(fvc::grad(alpha1_))//fluxScheme_->interpolate(alpha1_, "alpha")))
// //         );
// //
// //         forAll(alpha1_, celli)
// //         {
// //             if (alpha1_[celli] < 0.99 && alpha1_[celli] > 0.01)
// //             {
// //                 alpha1_[celli] +=
// //                     (2.0*alpha1_[celli] - 1.0)
// //                 *(1.0/Gamma[celli] - magGradAlpha[celli]);
// //             }
// //         }
// //         alpha1_.maxMin(0.0, 1.0);
//
//         scalar fr = 0.5;
//         alpha1_ = alpha1_*fr + (1.0 - fr)*sharpen(alpha1_);//H(levelSet_);
//     }
    alpha1_.correctBoundaryConditions();
    alpha2_ = 1.0 - alpha1_;




    // Update thermo
    thermo_.solve();

    // Surface tension force
    volVectorField sTF
    (
        fvc::reconstruct
        (
            interfaceProperties_.surfaceTensionForce()
           *rho_.mesh().magSf()
        )
    );

    //- Calculate deltas for momentum and energy
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
            fvm::ddt(alpha1_, rho1_) - fvc::ddt(alpha1_, rho1_)
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
            fvm::ddt(alpha2_, rho2_) - fvc::ddt(alpha2_, rho2_)
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
    rho_ = alphaRho1_ + alphaRho2_;

    compressibleBlastSystem::postUpdate();
}


void Foam::interfaceTwoPhaseCompressibleSystem::decode()
{
    // Calculate densities
    alpha1_.maxMin(0.0, 1.0);
    alpha1_.correctBoundaryConditions();
    alpha2_ = 1.0 - alpha1_;

    volScalarField alpha1
    (
        max(alpha1_, thermo_.thermo(0).residualAlpha())
    );
    volScalarField alpha2
    (
        max(alpha2_, thermo_.thermo(1).residualAlpha())
    );

    alphaRho1_.max(0);
    rho1_.ref() = alphaRho1_()/alpha1();
    rho1_.max(small);
    rho1_.correctBoundaryConditions();
    alphaRho1_.boundaryFieldRef() = alpha1_.boundaryField()*rho1_.boundaryField();

    alphaRho2_.max(0);
    rho2_.ref() = alphaRho2_()/alpha2();
    rho2_.max(small);
    rho2_.correctBoundaryConditions();
    alphaRho2_.boundaryFieldRef() = alpha2_.boundaryField()*rho2_.boundaryField();

    rho_ = alphaRho1_ + alphaRho2_;

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
