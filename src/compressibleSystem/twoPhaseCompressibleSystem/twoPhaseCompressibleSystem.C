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
    eos_(rho_, e_, p_, dict),
    alpha_(eos_.alpha()),
    rho1_(eos_.rho1()),
    rho2_(eos_.rho2()),
    alphaRho1_(eos_.alphaRho1()),
    alphaRho2_(eos_.alphaRho2()),
    alphaPhi_(eos_.alphaPhi()),
    alphaRhoPhi1_(eos_.alphaRhoPhi1()),
    alphaRhoPhi2_(eos_.alphaRhoPhi2())
{
    alphaRhoPhi1_ = fvc::interpolate(alphaRho1_)*phi_;
    alphaRhoPhi2_ = fvc::interpolate(alphaRho2_)*phi_;
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
    volScalarField alphaOld(ai[stepi - 1]*alpha_);
    volScalarField alphaRho1Old(ai[stepi - 1]*alphaRho1_);
    volScalarField alphaRho2Old(ai[stepi - 1]*alphaRho2_);
    if (oldIs_[stepi - 1] != -1)
    {
        alphaOld_[oldIs_[stepi - 1]] = alpha_;
        alphaRho1Old_[oldIs_[stepi - 1]] = alphaRho1_;
        alphaRho2Old_[oldIs_[stepi - 1]] = alphaRho2_;
    }

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

    volScalarField deltaAlpha(fvc::div(alphaPhi_) - alpha_*fvc::div(phi_));
    volScalarField deltaAlphaRho1(fvc::div(alphaRhoPhi1_));
    volScalarField deltaAlphaRho2(fvc::div(alphaRhoPhi2_));
    if (deltaIs_[stepi - 1] != -1)
    {
        deltaAlpha_[deltaIs_[stepi - 1]] = deltaAlpha;
        deltaAlphaRho1_[deltaIs_[stepi - 1]] = deltaAlphaRho1;
        deltaAlphaRho2_[deltaIs_[stepi - 1]] = deltaAlphaRho2;
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
    alpha_ = alphaOld - dT*deltaAlpha;
    alpha_.correctBoundaryConditions();

    alphaRho1_ = alphaRho1Old - dT*deltaAlphaRho1;
    alphaRho1_.correctBoundaryConditions();

    alphaRho2_ = alphaRho2Old - dT*deltaAlphaRho2;
    alphaRho2_.correctBoundaryConditions();

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
    for (label i = 0; i < nSteps; i++)
    {
        if (storeFields[i])
        {
            alphaOld_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName(alpha_.name(), Foam::name(i)),
                        rho_.time().timeName(),
                        rho_.mesh()
                    ),
                    alpha_
                )
            );
            alphaRho1Old_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName(alphaRho1_.name(), Foam::name(i)),
                        rho_.time().timeName(),
                        rho_.mesh()
                    ),
                    alphaRho1_
                )
            );
            alphaRho2Old_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName(alphaRho2_.name(), Foam::name(i)),
                        rho_.time().timeName(),
                        rho_.mesh()
                    ),
                    alphaRho2_
                )
            );
        }
    }
    for (label i = 0; i <+ nSteps; i++)
    {
        if (storeDeltas[i])
        {
            deltaAlpha_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            alpha_.name() + "Delta", Foam::name(i)
                        ),
                        rho_.time().timeName(),
                        rho_.mesh()
                    ),
                    rho_.mesh(),
                    dimensionedScalar("0", inv(dimTime), 0.0)
                )
            );
            deltaAlphaRho1_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            alphaRho1_.name() + "Delta", Foam::name(i)
                        ),
                        rho_.time().timeName(),
                        rho_.mesh()
                    ),
                    rho_.mesh(),
                    dimensionedScalar("0", dimDensity/dimTime, 0.0)
                )
            );
            deltaAlphaRho2_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            alphaRho2_.name() + "Delta", Foam::name(i)
                        ),
                        rho_.time().timeName(),
                        rho_.mesh()
                    ),
                    rho_.mesh(),
                    dimensionedScalar("0", dimDensity/dimTime, 0.0)
                )
            );
        }
    }
}


void Foam::twoPhaseCompressibleSystem::update()
{
    fluxScheme_->update
    (
        alpha_,
        rho1_,
        rho2_,
        U_,
        e_,
        p_,
        eos_.c()(),
        phi_,
        alphaPhi_,
        alphaRhoPhi1_,
        alphaRhoPhi2_,
        rhoPhi_,
        rhoUPhi_,
        rhoEPhi_
    );
}

void Foam::twoPhaseCompressibleSystem::decode()
{
    alpha_.max(0);
    alpha_.min(1);
    alpha_.correctBoundaryConditions();

    volScalarField alpha1(max(alpha_, 1e-10));
    volScalarField alpha2(max(1.0 - alpha1, 1e-10));

    alphaRho1_.max(0);
    rho1_.ref() = alphaRho1_()/alpha1();
    rho1_.correctBoundaryConditions();
    alphaRho1_.boundaryFieldRef() =
        alpha_.boundaryField()*rho1_.boundaryField();

    alphaRho2_.max(0);
    rho2_.ref() = alphaRho2_()/alpha2();
    rho2_.correctBoundaryConditions();
    alphaRho2_.boundaryFieldRef() =
        (1.0 - alpha_.boundaryField())*rho2_.boundaryField();

    rho_ = alphaRho1_ + alphaRho2_;

    U_.ref() = rhoU_()/rho_();
    U_.correctBoundaryConditions();

    rhoU_.boundaryFieldRef() = rho_.boundaryField()*U_.boundaryField();

    volScalarField E(rhoE_/rho_);
    e_.ref() = E() - 0.5*magSqr(U_());

    //--- Hard limit, e
    if(min(e_).value() < 0)
    {
        WarningInFunction<< "Limiting e, min(e) = " << min(e_).value() << endl;
        e_.max(small);
        rhoE_.ref() = rho_()*(e_() + 0.5*magSqr(U_()));
    }
    e_.correctBoundaryConditions();

    rhoE_.boundaryFieldRef() =
        rho_.boundaryField()
       *(
            e_.boundaryField()
          + 0.5*magSqr(U_.boundaryField())
        );

    eos_.updateP();
}


void Foam::twoPhaseCompressibleSystem::encode()
{
    alphaRho1_ = alpha_*rho1_;
    alphaRho2_ = (1.0 - alpha_)*rho2_;
    rho_ = alphaRho1_ + alphaRho2_;
    rhoU_ = rho_*U_;
    rhoE_ = rho_*(e_ + 0.5*magSqr(U_));
}

// ************************************************************************* //
