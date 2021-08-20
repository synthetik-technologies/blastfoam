/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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

#include "twoPhaseCompressibleSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(twoPhaseCompressibleSystem, 0);
    addToRunTimeSelectionTable
    (
        phaseCompressibleSystem,
        twoPhaseCompressibleSystem,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseCompressibleSystem::twoPhaseCompressibleSystem
(
    const fvMesh& mesh
)
:
    phaseCompressibleSystem(2, mesh),
    thermo_
    (
        refCast<twoPhaseFluidBlastThermo>(thermoPtr_())
    ),
    volumeFraction_(thermo_.volumeFraction()),
    rho1_(thermo_.thermo(0).rho()),
    rho2_(thermo_.thermo(1).rho()),
    alphaRho1_
    (
        IOobject
        (
            IOobject::groupName("alphaRho", rho1_.group()),
            mesh.time().timeName(),
            mesh
        ),
        volumeFraction_*rho1_
    ),
    alphaRho2_
    (
        IOobject
        (
            IOobject::groupName("alphaRho", rho2_.group()),
            mesh.time().timeName(),
            mesh
        ),
        (1.0 - volumeFraction_)*rho2_
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
    )
{
    rho_ = alphaRho1_ + alphaRho2_;

    thermo_.initializeModels();
    this->setModels();

    encode();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseCompressibleSystem::~twoPhaseCompressibleSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::twoPhaseCompressibleSystem::solve()
{
    //- Update changes in volume fraction and phase mass
    volScalarField deltaAlpha
    (
        fvc::div(alphaPhi_)
      - volumeFraction_*fvc::div(phi_)
    );
    volScalarField deltaAlphaRho1(fvc::div(alphaRhoPhi1_));
    volScalarField deltaAlphaRho2(fvc::div(alphaRhoPhi2_));

    this->storeAndBlendDelta(deltaAlpha);
    this->storeAndBlendDelta(deltaAlphaRho1);
    this->storeAndBlendDelta(deltaAlphaRho2);


    dimensionedScalar dT = rho_.time().deltaT();

    // Volume fraction is not scaled by change in volume because it is not
    // conserved
    this->storeAndBlendOld(volumeFraction_, false);
    this->storeAndBlendOld(alphaRho1_);
    this->storeAndBlendOld(alphaRho2_);
    rho_ = alphaRho1_ + alphaRho2_;
    rho_.storePrevIter();

    volumeFraction_ -= dT*deltaAlpha;
    volumeFraction_.correctBoundaryConditions();

    alphaRho1_.storePrevIter();
    alphaRho1_ -= dT*deltaAlphaRho1;
    alphaRho1_.correctBoundaryConditions();

    alphaRho2_.storePrevIter();
    alphaRho2_ -= dT*deltaAlphaRho2;
    alphaRho2_.correctBoundaryConditions();

    rho_ = alphaRho1_ + alphaRho2_;

    thermo_.solve();
    phaseCompressibleSystem::solve();
}


void Foam::twoPhaseCompressibleSystem::update()
{
    decode();
    fluxScheme_->update
    (
        volumeFraction_,
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
    thermo_.update();
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseCompressibleSystem::ESource() const
{
    return thermo_.ESource();
}

void Foam::twoPhaseCompressibleSystem::calcAlphaAndRho()
{
    volumeFraction_.max(0);
    volumeFraction_.min(1);
    volumeFraction_.correctBoundaryConditions();

    volScalarField alpha1(max(volumeFraction_, thermo_.thermo(0).residualAlpha()));
    volScalarField alpha2(max(1.0 - volumeFraction_, thermo_.thermo(1).residualAlpha()));

    alphaRho1_.max(0);
    rho1_.ref() = alphaRho1_()/alpha1();
    rho1_.max(small);
    rho1_.correctBoundaryConditions();
    alphaRho1_ = volumeFraction_*rho1_;

    alphaRho2_.max(0);
    rho2_.ref() = alphaRho2_()/alpha2();
    rho2_.max(small);
    rho2_.correctBoundaryConditions();
    alphaRho2_ = (1.0 - volumeFraction_)*rho2_;

    rho_ = alphaRho1_ + alphaRho2_;
}


void Foam::twoPhaseCompressibleSystem::decode()
{
    calcAlphaAndRho();

    U_.ref() = rhoU_()/rho_();
    U_.correctBoundaryConditions();

    rhoU_.boundaryFieldRef() = rho_.boundaryField()*U_.boundaryField();

    volScalarField E(rhoE_/rho_);
    e_.ref() = E() - 0.5*magSqr(U_());

    //- Limit internal energy it there is a negative temperature
    if(min(T_).value() < TLow_.value() && thermo_.limit())
    {
        if (debug)
        {
            WarningInFunction
                << "Lower limit of temperature reached, min(T) = "
                << min(T_).value()
                << ", limiting internal energy." << endl;
        }
        volScalarField limit(pos(T_ - TLow_));
        T_.max(TLow_);
        e_ = e_*limit + thermo_.he(p_, T_)*(1.0 - limit);
        rhoE_.ref() = rho_*(e_() + 0.5*magSqr(U_()));
    }
    e_.correctBoundaryConditions();

    rhoE_.boundaryFieldRef() =
        rho_.boundaryField()
       *(
            e_.boundaryField()
          + 0.5*magSqr(U_.boundaryField())
        );

    thermo_.correct();
    if (constraints_.constrainsField(p_.name()))
    {
        constraints_.constrain(p_);
        p_.correctBoundaryConditions();
    }
}


void Foam::twoPhaseCompressibleSystem::encode()
{
    alphaRho1_ = volumeFraction_*rho1_;
    alphaRho2_ = (1.0 - volumeFraction_)*rho2_;
    rho_ = alphaRho1_ + alphaRho2_;
    rhoU_ = rho_*U_;
    rhoE_ = rho_*(e_ + 0.5*magSqr(U_));
}

// ************************************************************************* //
