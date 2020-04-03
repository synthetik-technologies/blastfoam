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

#include "reactingCompressibleSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(reactingCompressibleSystem, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactingCompressibleSystem::reactingCompressibleSystem
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    phaseCompressibleSystem(mesh, dict),
    thermo_(rhoReactionThermo::New(mesh)),
    e_(thermo_->he())
{
    thermo_->validate("compressibleSystem", "e");
    rho_ = thermo_->rho();

    setModels(dict);

    Switch useChemistry
    (
        word(thermo_->subDict("thermoType").lookup("mixture"))
     == "reactingMixture"
    );
    if (useChemistry)
    {
        reaction_.set
        (
            CombustionModel<rhoReactionThermo>::New
            (
                refCast<rhoReactionThermo>(thermo_()),
                turbulence_()
            ).ptr()
        );
        word inertSpecie(thermo_->lookup("inertSpecie"));
        inertIndex_ = thermo_->composition().species()[inertSpecie];
    }
    encode();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reactingCompressibleSystem::~reactingCompressibleSystem()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::reactingCompressibleSystem::solve
(
    const label stepi,
    const scalarList& ai,
    const scalarList& bi
)
{
    if (oldIs_[stepi - 1] != -1)
    {
        rhoOld_.set
        (
            oldIs_[stepi - 1],
            new volScalarField(rho_)
        );
    }

    volScalarField rhoOld(ai[stepi - 1]*rho_);
    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = oldIs_[i];
        if (fi != -1 && ai[fi] != 0)
        {
            rhoOld += ai[fi]*rhoOld_[fi];
        }
    }

    volScalarField deltaRho(fvc::div(rhoPhi_));
    if (deltaIs_[stepi - 1] != -1)
    {
        deltaRho_.set
        (
            deltaIs_[stepi - 1],
            new volScalarField(deltaRho)
        );
    }
    deltaRho *= bi[stepi - 1];

    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = deltaIs_[i];
        if (fi != -1 && bi[fi] != 0)
        {
            deltaRho += bi[fi]*deltaRho_[fi];
        }
    }

    dimensionedScalar dT = rho_.time().deltaT();
    rho_ = rhoOld - dT*deltaRho;
    rho_.correctBoundaryConditions();

    phaseCompressibleSystem::solve(stepi, ai, bi);

    //- Solve chemistry on the final step
    if (stepi == oldIs_.size())
    {
        if (reaction_.valid())
        {
            reaction_->correct();

            PtrList<volScalarField>& Y = thermo_->composition().Y();
            volScalarField Yt(0.0*Y[0]);
            forAll(Y, i)
            {
                if (i != inertIndex_ && thermo_->composition().active(i))
                {
                    volScalarField& Yi = Y[i];

                    fvScalarMatrix YiEqn
                    (
                        fvm::ddt(rho_, Yi)
                      + fvm::div(rhoPhi_, Yi, "div(rhoPhi,Yi)")
                      - fvm::laplacian(turbulence_->alphaEff(), Yi)
                    ==
                        reaction_->R(Yi)
                    );
                    YiEqn.solve("Yi");


                    Yi.max(0.0);
                    Yt += Yi;
                }
            }
            Y[inertIndex_] = scalar(1) - Yt;
            Y[inertIndex_].max(0.0);
            rhoE_ += dT*reaction_->Qdot();
        }
    }
    decode();
}


void Foam::reactingCompressibleSystem::setODEFields
(
    const label nSteps,
    const boolList& storeFields,
    const boolList& storeDeltas
)
{
    phaseCompressibleSystem::setODEFields(nSteps, storeFields, storeDeltas);
    rhoOld_.setSize(nOld_);

    deltaRho_.setSize(nDelta_);
}

void Foam::reactingCompressibleSystem::clearODEFields()
{
    phaseCompressibleSystem::clearODEFields();

    rhoOld_.clear();
    rhoOld_.setSize(nOld_);

    deltaRho_.clear();
    deltaRho_.setSize(nDelta_);
}


void Foam::reactingCompressibleSystem::update()
{
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
}


Foam::tmp<Foam::volScalarField>
Foam::reactingCompressibleSystem::ESource() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "ESource",
                rho_.mesh().time().timeName(),
                rho_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            rho_.mesh(),
            dimensionedScalar("0", rhoE_.dimensions()/dimTime, 0.0)
        )
    );
}


void Foam::reactingCompressibleSystem::decode()
{
    U_.ref() = rhoU_()/rho_();
    U_.correctBoundaryConditions();

    rhoU_.boundaryFieldRef() = rho_.boundaryField()*U_.boundaryField();

    volScalarField E(rhoE_/rho_);
    e_.ref() = E() - 0.5*magSqr(U_());
    e_.correctBoundaryConditions();

    rhoE_.boundaryFieldRef() =
        rho_.boundaryField()
       *(
            e_.boundaryField()
          + 0.5*magSqr(U_.boundaryField())
        );

    thermo_->correct();
    p_.ref() = rho_/thermo_->psi();
    p_.correctBoundaryConditions();
    rho_.boundaryFieldRef() ==
        thermo_->psi().boundaryField()*p_.boundaryField();
}


void Foam::reactingCompressibleSystem::encode()
{
    rhoU_ = rho_*U_;
    rhoE_ = rho_*(e_ + 0.5*magSqr(U_));
}


Foam::tmp<Foam::volScalarField>
Foam::reactingCompressibleSystem::speedOfSound() const
{
    return sqrt(thermo_->gamma()/thermo_->psi());
}


Foam::tmp<Foam::volScalarField> Foam::reactingCompressibleSystem::Cv() const
{
    return thermo_->Cv();
}


Foam::tmp<Foam::volScalarField> Foam::reactingCompressibleSystem::mu() const
{
    return thermo_->mu();
}


Foam::tmp<Foam::scalarField>
Foam::reactingCompressibleSystem::mu(const label patchi) const
{
    return thermo_->mu(patchi);
}

Foam::tmp<Foam::volScalarField> Foam::reactingCompressibleSystem::nu() const
{
    return thermo_->nu();
}

Foam::tmp<Foam::scalarField>
Foam::reactingCompressibleSystem::nu(const label patchi) const
{
    return thermo_->nu(patchi);
}

Foam::tmp<Foam::volScalarField>
Foam::reactingCompressibleSystem::alpha() const
{
    return thermo_->alpha();
}

Foam::tmp<Foam::scalarField>
Foam::reactingCompressibleSystem::alpha(const label patchi) const
{
    return thermo_->alpha(patchi);
}

Foam::tmp<Foam::volScalarField> Foam::reactingCompressibleSystem::alphaEff
(
    const volScalarField& alphat
) const
{
    return thermo_->alphaEff(alphat);
}

Foam::tmp<Foam::scalarField> Foam::reactingCompressibleSystem::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return thermo_->alphaEff(alphat, patchi);
}

Foam::tmp<Foam::volScalarField>
Foam::reactingCompressibleSystem::alphahe() const
{
    return thermo_->alphahe();
}

Foam::tmp<Foam::scalarField>
Foam::reactingCompressibleSystem::alphahe(const label patchi) const
{
    return thermo_->alphahe(patchi);
}

Foam::tmp<Foam::volScalarField> Foam::reactingCompressibleSystem::kappa() const
{
    return thermo_->kappa();
}

Foam::tmp<Foam::scalarField>
Foam::reactingCompressibleSystem::kappa(const label patchi) const
{
    return thermo_->kappa(patchi);
}

Foam::tmp<Foam::volScalarField> Foam::reactingCompressibleSystem::kappaEff
(
    const volScalarField& alphat
) const
{
    return thermo_->kappaEff(alphat);
}

Foam::tmp<Foam::scalarField> Foam::reactingCompressibleSystem::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return thermo_->kappaEff(alphat, patchi);
}
// ************************************************************************* //
