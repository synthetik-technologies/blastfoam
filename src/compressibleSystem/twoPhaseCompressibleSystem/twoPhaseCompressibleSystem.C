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
    const fvMesh& mesh,
    const dictionary& dict
)
:
    phaseCompressibleSystem(mesh, dict),
    thermo_(word::null, p_, rho_, e_, T_, dict, true),
    volumeFraction_(thermo_.volumeFraction()),
    rho1_(thermo_.thermo1().rho()),
    rho2_(thermo_.thermo2().rho()),
    alphaRho1_
    (
        IOobject
        (
            IOobject::groupName("alphaRho", rho1_.group()),
            mesh.time().timeName(),
            mesh
        ),
        volumeFraction_*rho1_,
        wordList(p_.boundaryField().types().size(), "zeroGradient")
    ),
    alphaRho2_
    (
        IOobject
        (
            IOobject::groupName("alphaRho", rho2_.group()),
            mesh.time().timeName(),
            mesh
        ),
        (1.0 - volumeFraction_)*rho2_,
        wordList(p_.boundaryField().types().size(), "zeroGradient")
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
    encode();
    setModels(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseCompressibleSystem::~twoPhaseCompressibleSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::twoPhaseCompressibleSystem::solve
(
    const label stepi,
    const scalarList& ai,
    const scalarList& bi
)
{
    if (oldIs_[stepi - 1] != -1)
    {
        alphaOld_.set
        (
            oldIs_[stepi - 1],
            new volScalarField(volumeFraction_)
        );
        alphaRho1Old_.set
        (
            oldIs_[stepi - 1],
            new volScalarField(alphaRho1_)
        );
        alphaRho2Old_.set
        (
            oldIs_[stepi - 1],
            new volScalarField(alphaRho2_)
        );
    }

    volScalarField alphaOld(ai[stepi - 1]*volumeFraction_);
    volScalarField alphaRho1Old(ai[stepi - 1]*alphaRho1_);
    volScalarField alphaRho2Old(ai[stepi - 1]*alphaRho2_);
    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = oldIs_[i];
        if (fi != -1 && ai[fi] != 0)
        {
            alphaOld += ai[fi]*alphaOld_[fi];
            alphaRho1Old += ai[fi]*alphaRho1Old_[fi];
            alphaRho2Old += ai[fi]*alphaRho2Old_[fi];
        }
    }

    volScalarField deltaAlpha
    (
        fvc::div(alphaPhi_)
      - volumeFraction_*fvc::div(phi_)
    );
    volScalarField deltaAlphaRho1(fvc::div(alphaRhoPhi1_));
    volScalarField deltaAlphaRho2(fvc::div(alphaRhoPhi2_));
    if (deltaIs_[stepi - 1] != -1)
    {
        deltaAlpha_.set
        (
            deltaIs_[stepi - 1],
            new volScalarField(deltaAlpha)
        );
        deltaAlphaRho1_.set
        (
            deltaIs_[stepi - 1],
            new volScalarField(deltaAlphaRho1)
        );
        deltaAlphaRho2_.set
        (
            deltaIs_[stepi - 1],
            new volScalarField(deltaAlphaRho2)
        );
    }
    deltaAlpha *= bi[stepi - 1];
    deltaAlphaRho1 *= bi[stepi - 1];
    deltaAlphaRho2 *= bi[stepi - 1];

    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = deltaIs_[i];
        if (fi != -1 && bi[fi] != 0)
        {
            deltaAlpha += bi[fi]*deltaAlpha_[fi];
            deltaAlphaRho1 += bi[fi]*deltaAlphaRho1_[fi];
            deltaAlphaRho2 += bi[fi]*deltaAlphaRho2_[fi];
        }
    }

    dimensionedScalar dT = rho_.time().deltaT();
    volumeFraction_ = alphaOld - dT*deltaAlpha;
    volumeFraction_.correctBoundaryConditions();

    alphaRho1_.oldTime() = alphaRho1Old;
    alphaRho1_ = alphaRho1Old - dT*deltaAlphaRho1;
    alphaRho1_.correctBoundaryConditions();

    alphaRho2_.oldTime() = alphaRho2Old;
    alphaRho2_ = alphaRho2Old - dT*deltaAlphaRho2;
    alphaRho2_.correctBoundaryConditions();

    thermo_.solve(stepi, ai, bi);
    phaseCompressibleSystem::solve(stepi, ai, bi);
    decode();
}


void Foam::twoPhaseCompressibleSystem::setODEFields
(
    const label nSteps,
    const boolList& storeFields,
    const boolList& storeDeltas
)
{
    phaseCompressibleSystem::setODEFields(nSteps, storeFields, storeDeltas);
    alphaOld_.setSize(nOld_);
    alphaRho1Old_.setSize(nOld_);
    alphaRho2Old_.setSize(nOld_);

    deltaAlpha_.setSize(nDelta_);
    deltaAlphaRho1_.setSize(nDelta_);
    deltaAlphaRho2_.setSize(nDelta_);

    thermo_.setODEFields(nSteps, oldIs_, nOld_, deltaIs_, nDelta_);
}

void Foam::twoPhaseCompressibleSystem::clearODEFields()
{
    phaseCompressibleSystem::clearODEFields();

    alphaOld_.clear();
    alphaRho1Old_.clear();
    alphaRho2Old_.clear();

    alphaOld_.setSize(nOld_);
    alphaRho1Old_.setSize(nOld_);
    alphaRho2Old_.setSize(nOld_);

    deltaAlpha_.clear();
    deltaAlphaRho1_.clear();
    deltaAlphaRho2_.clear();

    deltaAlpha_.setSize(nDelta_);
    deltaAlphaRho1_.setSize(nDelta_);
    deltaAlphaRho2_.setSize(nDelta_);

    thermo_.clearODEFields();
}


void Foam::twoPhaseCompressibleSystem::update()
{
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

    volScalarField alpha1(max(volumeFraction_, thermo_.thermo1().residualAlpha()));
    volScalarField alpha2(max(1.0 - volumeFraction_, thermo_.thermo2().residualAlpha()));

    alphaRho1_.max(0);
    rho1_.ref() = alphaRho1_()/alpha1();
    rho1_.correctBoundaryConditions();
    alphaRho1_ = volumeFraction_*rho1_;

    alphaRho2_.max(0);
    rho2_.ref() = alphaRho2_()/alpha2();
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
        e_ = e_*limit + thermo_.E()*(1.0 - limit);
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
}


void Foam::twoPhaseCompressibleSystem::encode()
{
    volumeFraction_.correctBoundaryConditions();
    rho1_.correctBoundaryConditions();
    rho2_.correctBoundaryConditions();
    U_.correctBoundaryConditions();
    e_.correctBoundaryConditions();

    alphaRho1_ = volumeFraction_*rho1_;
    alphaRho2_ = (1.0 - volumeFraction_)*rho2_;
    rho_ = alphaRho1_ + alphaRho2_;
    rhoU_ = rho_*U_;
    rhoE_ = rho_*(e_ + 0.5*magSqr(U_));
}

// ************************************************************************* //
