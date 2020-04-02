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
    eos_(rho_, e_, p_, dict),
    alphas_(eos_.alphas()),
    rhos_(eos_.rhos()),
    alphaRhos_(eos_.alphaRhos()),
    alphaPhis_(eos_.alphaPhis()),
    alphaRhoPhis_(eos_.alphaRhoPhis())
{
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
            phasei, new volScalarField(ai[stepi - 1]*alphas_[phasei])
        );
        alphaRhosOld.set
        (
            phasei, new volScalarField(ai[stepi - 1]*alphaRhos_[phasei])
        );
    }
    if (oldIs_[stepi - 1] != -1)
    {
        forAll(alphas_, phasei)
        {
            alphasOld_[oldIs_[stepi - 1]][phasei] = alphas_[phasei];
            alphaRhosOld_[oldIs_[stepi - 1]][phasei] = alphaRhos_[phasei];
        }
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
            deltaAlphas_[deltaIs_[stepi - 1]][phasei] = deltaAlphas[phasei];
            deltaAlphaRhos_[deltaIs_[stepi - 1]][phasei] =
                deltaAlphaRhos[phasei];
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

        alphaRhos_[phasei] = alphaRhosOld[phasei] - dT*deltaAlphaRhos[phasei];
        alphaRhos_[phasei].correctBoundaryConditions();
    }
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
    label nOlds = 0;
    forAll(oldIs_, i)
    {
        if (oldIs_[i])
        {
            nOlds++;
        }
    }
    alphasOld_.resize(nOlds);
    alphaRhosOld_.resize(nOlds);
    label step = 0;
    for (label i = 0; i < nSteps; i++)
    {
        if (storeFields[i])
        {
            alphasOld_.set(step, new PtrList<volScalarField>(alphas_.size()));
            alphaRhosOld_.set(step, new PtrList<volScalarField>(alphas_.size()));

            forAll(alphas_, phasei)
            {
                alphasOld_[step].set
                (
                    phasei,
                    new volScalarField
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                alphas_[phasei].name(), Foam::name(i)
                            ),
                            rho_.time().timeName(),
                            rho_.mesh()
                        ),
                        alphas_[phasei]
                    )
                );
                alphaRhosOld_[step].set
                (
                    phasei,
                    new volScalarField
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                alphaRhos_[phasei].name(), Foam::name(i)
                            ),
                            rho_.time().timeName(),
                            rho_.mesh()
                        ),
                        alphaRhos_[phasei]
                    )
                );
            }
            step++;
        }
    }

    label nDeltas = 0;
    forAll(deltaIs_, i)
    {
        if (deltaIs_[i])
        {
            nDeltas++;
        }
    }
    deltaAlphas_.resize(nDeltas);
    deltaAlphaRhos_.resize(nDeltas);
    step = 0;
    for (label i = 0; i <+ nSteps; i++)
    {
        if (storeDeltas[i])
        {
            deltaAlphas_.set(step, new PtrList<volScalarField>(alphas_.size()));
            deltaAlphaRhos_.set
            (
                step, new PtrList<volScalarField>(alphas_.size())
            );
            forAll(alphas_, phasei)
            {
                deltaAlphas_[step].set
                (
                    phasei,
                    new volScalarField
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                alphas_[phasei].name() + "Delta", Foam::name(i)
                            ),
                            rho_.time().timeName(),
                            rho_.mesh()
                        ),
                        rho_.mesh(),
                        dimensionedScalar("0", inv(dimTime), 0.0)
                    )
                );
                deltaAlphaRhos_[step].set
                (
                    phasei,
                    new volScalarField
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                alphaRhos_[phasei].name() + "Delta",
                                Foam::name(i)
                            ),
                            rho_.time().timeName(),
                            rho_.mesh()
                        ),
                        rho_.mesh(),
                        dimensionedScalar("0", dimDensity/dimTime, 0.0)
                    )
                );
            }
            step++;
        }
    }
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
        eos_.c()(),
        phi_,
        alphaPhis_,
        alphaRhoPhis_,
        rhoPhi_,
        rhoUPhi_,
        rhoEPhi_
    );
}

void Foam::multiphaseCompressibleSystem::decode()
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
    forAll(alphas_, phasei)
    {
        alphas_[phasei].max(0);
        alphas_[phasei].min(1);
        alphas_[phasei].correctBoundaryConditions();
        sumAlpha += alphas_[phasei];

        alphaRhos_[phasei].max(0);
        rhos_[phasei] = alphaRhos_[phasei]/max(alphas_[phasei], 1e-6);
        rhos_[phasei].correctBoundaryConditions();
    }
    Info<<"sumAlpha max = " << max(sumAlpha).value()
        <<", min = " << min(sumAlpha).value() << endl;

    forAll(alphas_, phasei)
    {
        alphas_[phasei] /= sumAlpha;
        alphas_[phasei].correctBoundaryConditions();

        alphaRhos_[phasei] = alphas_[phasei]*rhos_[phasei];
        rho_ += alphaRhos_[phasei];
    }

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
