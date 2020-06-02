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

#include "multiphaseCompressibleSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiphaseCompressibleSystem, 0);
    addToRunTimeSelectionTable
    (
        phaseCompressibleSystem,
        multiphaseCompressibleSystem,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseCompressibleSystem::multiphaseCompressibleSystem
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    phaseCompressibleSystem(mesh, dict),
    thermo_(word::null, p_, rho_, e_, T_, dict, true),
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

    setModels(dict);
    encode();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiphaseCompressibleSystem::~multiphaseCompressibleSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multiphaseCompressibleSystem::solve
(
    const label stepi,
    const scalarList& ai,
    const scalarList& bi
)
{
    PtrList<volScalarField> alphasOld(alphas_.size());
    PtrList<volScalarField> alphaRhosOld(alphas_.size());
    forAll(alphas_, phasei)
    {
        alphasOld.set
        (
            phasei, new volScalarField(alphas_[phasei])
        );
        alphaRhosOld.set
        (
            phasei, new volScalarField(alphaRhos_[phasei])
        );
    }
    if (oldIs_[stepi - 1] != -1)
    {
        forAll(alphas_, phasei)
        {
            alphasOld_[oldIs_[stepi - 1]].set
            (
                phasei,
                new volScalarField(alphas_[phasei])
            );
            alphaRhosOld_[oldIs_[stepi - 1]].set
            (
                phasei,
                new volScalarField(alphaRhos_[phasei])
            );
        }
    }

    forAll(alphas_, phasei)
    {
        alphasOld[phasei] *= ai[stepi - 1];
        alphaRhosOld[phasei] *= ai[stepi - 1];
    }
    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = oldIs_[i];
        if (fi != -1 && ai[fi] != 0)
        {
            forAll(alphas_, phasei)
            {
                alphasOld[phasei] += ai[fi]*alphasOld_[fi][phasei];
                alphaRhosOld[phasei] += ai[fi]*alphaRhosOld_[fi][phasei];
            }
        }
    }

    PtrList<volScalarField> deltaAlphas(alphas_.size());
    PtrList<volScalarField> deltaAlphaRhos(alphas_.size());
    forAll(alphas_, phasei)
    {
        deltaAlphas.set
        (
            phasei,
            new volScalarField
            (
                fvc::div(alphaPhis_[phasei])
              - alphas_[phasei]*fvc::div(phi_)
            )
        );
        deltaAlphaRhos.set
        (
            phasei, new volScalarField(fvc::div(alphaRhoPhis_[phasei]))
        );
    }
    if (deltaIs_[stepi - 1] != -1)
    {
        forAll(alphas_, phasei)
        {
            deltaAlphas_[deltaIs_[stepi - 1]].set
            (
                phasei,
                new volScalarField(deltaAlphas[phasei])
            );
            deltaAlphaRhos_[deltaIs_[stepi - 1]].set
            (
                phasei,
                new volScalarField(deltaAlphaRhos[phasei])
            );
        }
    }
    forAll(alphas_, phasei)
    {
        deltaAlphas[phasei] *= bi[stepi - 1];
        deltaAlphaRhos[phasei] *= bi[stepi - 1];
    }

    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = deltaIs_[i];
        if (fi != -1 && bi[fi] != 0)
        {
            forAll(alphas_, phasei)
            {
                deltaAlphas[phasei] += bi[fi]*deltaAlphas_[fi][phasei];
                deltaAlphaRhos[phasei] += bi[fi]*deltaAlphaRhos_[fi][phasei];
            }
        }
    }

    dimensionedScalar dT = rho_.time().deltaT();

    forAll(alphas_, phasei)
    {
        alphas_[phasei] = alphasOld[phasei] - dT*deltaAlphas[phasei];
        alphas_[phasei].correctBoundaryConditions();

        alphaRhos_[phasei].oldTime() = alphaRhosOld[phasei];
        alphaRhos_[phasei] = alphaRhosOld[phasei] - dT*deltaAlphaRhos[phasei];
        alphaRhos_[phasei].correctBoundaryConditions();
    }

    thermo_.solve(stepi, ai, bi);
    phaseCompressibleSystem::solve(stepi, ai, bi);

    decode();
}


void Foam::multiphaseCompressibleSystem::setODEFields
(
    const label nSteps,
    const boolList& storeFields,
    const boolList& storeDeltas
)
{
    phaseCompressibleSystem::setODEFields(nSteps, storeFields, storeDeltas);

    alphasOld_.resize(nOld_);
    alphaRhosOld_.resize(nOld_);
    forAll(alphasOld_, stepi)
    {
        alphasOld_.set
        (
            stepi,
            new PtrList<volScalarField>(alphas_.size())
        );
        alphaRhosOld_.set
        (
            stepi,
            new PtrList<volScalarField>(alphas_.size())
        );
    }

    deltaAlphas_.resize(nDelta_);
    deltaAlphaRhos_.resize(nDelta_);
    forAll(deltaAlphas_, stepi)
    {
        deltaAlphas_.set
        (
            stepi,
            new PtrList<volScalarField>(alphas_.size())
        );
        deltaAlphaRhos_.set
        (
            stepi,
            new PtrList<volScalarField>(alphas_.size())
        );
    }
    thermo_.setODEFields(nSteps, oldIs_, nOld_, deltaIs_, nDelta_);
}

void Foam::multiphaseCompressibleSystem::clearODEFields()
{
    phaseCompressibleSystem::clearODEFields();

    forAll(alphasOld_, stepi)
    {
        alphasOld_[stepi].clear();
        alphaRhosOld_[stepi].clear();
        alphasOld_[stepi].resize(alphas_.size());
        alphaRhosOld_[stepi].resize(alphas_.size());
    }

    forAll(deltaAlphas_, stepi)
    {
        deltaAlphas_[stepi].clear();
        deltaAlphaRhos_[stepi].clear();
        deltaAlphas_[stepi].resize(alphas_.size());
        deltaAlphaRhos_[stepi].resize(alphas_.size());
    }
    thermo_.clearODEFields();
}

void Foam::multiphaseCompressibleSystem::update()
{
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
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseCompressibleSystem::ESource() const
{
    return thermo_.ESource();
}


void Foam::multiphaseCompressibleSystem::calcAlphaAndRho()
{
    rho_ = dimensionedScalar("0", dimDensity, 0.0);
    volScalarField sumAlpha
    (
        IOobject
        (
            "sumAlpha",
            rho_.time().timeName(),
            rho_.mesh()
        ),
        rho_.mesh(),
        0.0
    );
    for (label phasei = 0; phasei < alphas_.size() - 1; phasei++)
    {
        alphas_[phasei].max(0);
        alphas_[phasei].min(1);
        alphas_[phasei].correctBoundaryConditions();
        sumAlpha += alphas_[phasei];

        alphaRhos_[phasei].max(0);
        rhos_[phasei] =
            alphaRhos_[phasei]
           /max(alphas_[phasei], thermo_.thermo(phasei).residualAlpha());

        rhos_[phasei].correctBoundaryConditions();

        alphaRhos_[phasei] = alphas_[phasei]*rhos_[phasei];

        rho_ += alphaRhos_[phasei];
    }
    label phasei = alphas_.size() - 1;
    alphas_[phasei] = 1.0 - sumAlpha;
    alphas_[phasei].max(0.0);
    alphas_[phasei].min(1.0);

    alphaRhos_[phasei].max(0.0);
    rhos_[phasei] = alphaRhos_[phasei]
           /max(alphas_[phasei], thermo_.thermo(phasei).residualAlpha());
    rhos_[phasei].correctBoundaryConditions();
    alphaRhos_[phasei] = alphas_[phasei]*rhos_[phasei];
    rho_ += alphaRhos_[phasei];
}

void Foam::multiphaseCompressibleSystem::decode()
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


void Foam::multiphaseCompressibleSystem::encode()
{
    rho_ = dimensionedScalar("0", dimDensity, 0.0);
    forAll(alphas_, phasei)
    {
        alphaRhos_[phasei] = alphas_[phasei]*rhos_[phasei];
        rho_ += alphaRhos_[phasei];
    }
    rhoU_ = rho_*U_;
    rhoE_ = rho_*(e_ + 0.5*magSqr(U_));
}

// ************************************************************************* //
