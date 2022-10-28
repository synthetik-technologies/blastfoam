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
    phi_ = fvc::flux(U_);
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
    surfaceScalarField alpha2f(1.0 - alpha1f);

    alphaPhi_ = alpha1f*phi_;
    alphaRhoPhi1_ = fluxScheme_->flux(rho1_, alpha1f, alpha1f, phi_);
    alphaRhoPhi2_ = rhoPhi_ - alphaRhoPhi1_;
    thermo_.update();

    interfacePtr_->correct();
}

// ************************************************************************* //
