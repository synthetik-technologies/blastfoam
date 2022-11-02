/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2022
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

#include "twoPhaseInterfaceCompressibleSystem.H"
#include "addToRunTimeSelectionTable.H"
#include "EulerDdtScheme.H"
#include "gaussConvectionScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(twoPhaseInterfaceCompressibleSystem, 0);
    addToRunTimeSelectionTable
    (
        compressibleSystem,
        twoPhaseInterfaceCompressibleSystem,
        twoPhase
    );
}


// * * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::twoPhaseInterfaceCompressibleSystem::rhoUSource() const
{
    return
        compressibleBlastSystem::rhoUSource()
      + fvc::reconstruct(interfacePtr_->surfaceTensionForce()*mesh().magSf());
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseInterfaceCompressibleSystem::rhoESource() const
{
    return
        compressibleBlastSystem::rhoESource()
      + (
            fvc::reconstruct(interfacePtr_->surfaceTensionForce()*mesh().magSf())
          & U()
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseInterfaceCompressibleSystem::twoPhaseInterfaceCompressibleSystem
(
    const fvMesh& mesh
)
:
    twoPhaseCompressibleSystem(mesh),
    interfacePtr_
    (
        new interfaceProperties
        (
            alpha1_,
            alpha2_,
            U_,
            *this
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseInterfaceCompressibleSystem::~twoPhaseInterfaceCompressibleSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::twoPhaseInterfaceCompressibleSystem::update()
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

    surfaceScalarField alpha1f
    (
        fvc::interpolate(alpha1_, phi_, "reconstruct(alpha)")
    );

    tmp<surfaceScalarField> talphaRho1Own, talphaRho1Nei;
    tmp<surfaceScalarField> talphaRho2Own, talphaRho2Nei;
    {
        autoPtr<ReconstructionScheme<scalar>> rho1Limiter
        (
            ReconstructionScheme<scalar>::New(rho1_, "rho", rho1_.group())
        );

        autoPtr<ReconstructionScheme<scalar>> rho2Limiter
        (
            ReconstructionScheme<scalar>::New(rho2_, "rho", rho2_.group())
        );

        talphaRho1Own = surfaceScalarField::New
        (
            rho1Limiter->ownName(alphaRho1_.name()),
            alpha1f*rho1Limiter->interpolateOwn()
        );
        talphaRho1Nei = surfaceScalarField::New
        (
            rho1Limiter->neiName(alphaRho1_.name()),
            alpha1f*rho1Limiter->interpolateNei()
        );

        surfaceScalarField alpha2f(1.0 - alpha1f);
        talphaRho2Own = surfaceScalarField::New
        (
            rho2Limiter->ownName(alphaRho2_.name()),
            alpha2f*rho2Limiter->interpolateOwn()
        );
        talphaRho2Nei = surfaceScalarField::New
        (
            rho2Limiter->ownName(alphaRho2_.name()),
            alpha2f*rho2Limiter->interpolateNei()
        );
        static bool cached = false;
        if (!cached)
        {
            mesh().addTemporaryObject(talphaRho1Own().name());
            mesh().addTemporaryObject(talphaRho1Nei().name());
            mesh().addTemporaryObject(talphaRho2Own().name());
            mesh().addTemporaryObject(talphaRho2Nei().name());
        }
    }
    alphaPhi_ = alpha1f*phi_;
    alphaRhoPhi1_ =// phi_*(pos0(phi_)*talphaRho1Own + neg(phi_)*talphaRho1Nei);
        fluxScheme_->flux(talphaRho1Own(), talphaRho1Nei(), phi_);
    alphaRhoPhi2_ =// rhoPhi_ - alphaRhoPhi1_;
        fluxScheme_->flux(talphaRho2Own(), talphaRho2Nei(), phi_);
    thermo_.update();

    interfacePtr_->correct();
}

// ************************************************************************* //
