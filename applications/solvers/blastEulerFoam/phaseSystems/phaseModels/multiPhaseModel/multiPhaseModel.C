/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2019-04-29 Jeff Heylmun:    Simplified model
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "multiPhaseModel.H"
#include "phaseSystem.H"
#include "fvMatrix.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"
#include "phaseCompressibleTurbulenceModel.H"
#include "addToRunTimeSelectionTable.H"
#include "SortableList.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiPhaseModel, 0);
    addToRunTimeSelectionTable
    (
        phaseModel,
        multiPhaseModel,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiPhaseModel::multiPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const label index
)
:
    phaseModel(fluid, phaseName, index),
    phases_(phaseDict_.lookup("phases")),
    p_
    (
        IOobject
        (
            IOobject::groupName("p", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fluid_.p(),
        fluid_.p().boundaryField()
    ),
    rho_
    (
        IOobject
        (
            IOobject::groupName("rho", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh()
    ),
    alphaPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaPhi", name_),
            this->mesh().time().timeName(),
            this->mesh()
        ),
        this->mesh(),
        dimensionedScalar("0", phi_.dimensions(), 0.0)
    ),
    thermo_(phaseName, p_, rho_, e_, T_, phaseDict_, true, phaseName),
    alphas_(thermo_.volumeFractions()),
    rhos_(thermo_.rhos()),
    alphaRhos_(alphas_.size()),
    alphaPhis_(alphas_.size()),
    alphaRhoPhis_(alphas_.size()),
    alphasOld_(alphas_.size()),
    alphaRhosOld_(alphas_.size()),
    deltaAlphas_(alphas_.size()),
    deltaAlphaRhos_(alphas_.size()),
    fluxScheme_(fluxScheme::New(fluid.mesh(), name_))
{

    //- Temporarily Store read density
    volScalarField sumAlpha
    (
        IOobject
        (
            "sumAlpha",
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        0.0
    );
    alphaRho_ = dimensionedScalar(dimDensity, 0.0);

    forAll(alphas_, phasei)
    {
        sumAlpha += alphas_[phasei];
        word phaseName = alphas_[phasei].group();
        alphaRhos_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("alphaRho", phaseName),
                    fluid.mesh().time().timeName(),
                    fluid.mesh()
                ),
                alphas_[phasei]*rhos_[phasei],
                rhos_[phasei].boundaryField().types()
            )
        );
        alphaRho_ += alphaRhos_[phasei];
        alphaPhis_.set
        (
            phasei,
            new surfaceScalarField
            (
                IOobject
                (
                    IOobject::groupName("alphaPhi", phaseName),
                    fluid.mesh().time().timeName(),
                    fluid.mesh()
                ),
                fluid.mesh(),
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
                    fluid.mesh().time().timeName(),
                    fluid.mesh()
                ),
                fluid.mesh(),
                dimensionedScalar("0", dimensionSet(1, 0, -1, 0, 0), 0.0)
            )
        );

        alphasOld_.set
        (
            phasei,
            new PtrList<volScalarField>()
        );
        alphaRhosOld_.set
        (
            phasei,
            new PtrList<volScalarField>()
        );

        deltaAlphas_.set
        (
            phasei,
            new PtrList<volScalarField>()
        );
        deltaAlphaRhos_.set
        (
            phasei,
            new PtrList<volScalarField>()
        );
    }
    // Reset density to correct value
    volScalarField& alpha = *this;
    alpha = sumAlpha;
    alpha.correctBoundaryConditions();

    rho_ = alphaRho_/Foam::max(sumAlpha, residualAlpha());
    thermo_.setTotalVolumeFractionPtr(*this);

    this->turbulence_ =
        phaseCompressibleTurbulenceModel::New
        (
            *this,
            rho_,
            U_,
            alphaRhoPhi_,
            phi_,
            *this
        );
    phaseModel::initializeModels();
    thermo_.initializeModels();
    correctThermo();

    encode();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiPhaseModel::~multiPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multiPhaseModel::solveAlpha(const bool s)
{
    phaseModel::solveAlpha(s);
}


void Foam::multiPhaseModel::solve()
{
    dimensionedScalar dT = rho_.time().deltaT();
    forAll(alphas_, phasei)
    {
        volScalarField deltaAlpha
        (
            fvc::div(alphaPhis_[phasei]) - alphas_[phasei]*fvc::div(phi_)
        );
        this->storeAndBlendDelta(deltaAlpha, deltaAlphas_[phasei]);

        volScalarField deltaAlphaRho(fvc::div(alphaRhoPhis_[phasei]));
        this->storeAndBlendDelta(deltaAlphaRho, deltaAlphaRhos_[phasei]);

        this->storeAndBlendOld(alphas_[phasei], alphasOld_[phasei]);
        alphas_[phasei] -= dT*deltaAlpha;
        alphas_[phasei].correctBoundaryConditions();

        this->storeAndBlendOld(alphaRhos_[phasei], alphaRhosOld_[phasei]);
        alphaRhos_[phasei].storePrevIter();
        alphaRhos_[phasei] -= dT*deltaAlphaRho;
        alphaRhos_[phasei].correctBoundaryConditions();
    }

    thermo_.solve();
    phaseModel::solve();
}


void Foam::multiPhaseModel::postUpdate()
{
    phaseModel::postUpdate();
    thermo_.postUpdate();
}


void Foam::multiPhaseModel::clearODEFields()
{
    phaseModel::clearODEFields();
    fluxScheme_->clear();

    forAll(alphasOld_, phasei)
    {
        this->clearOld(alphasOld_[phasei]);
        this->clearOld(alphaRhosOld_[phasei]);

        this->clearDelta(deltaAlphas_[phasei]);
        this->clearDelta(deltaAlphaRhos_[phasei]);
    }
    thermo_.clearODEFields();
}


void Foam::multiPhaseModel::update()
{
    tmp<volScalarField> c(speedOfSound());
    fluxScheme_->update
    (
        alphas_,
        rhos_,
        U_,
        e_,
        p_,
        c(),
        phi_,
        alphaPhi_,
        alphaRhoPhi_,
        alphaPhis_,
        alphaRhoPhis_,
        alphaRhoUPhi_,
        alphaRhoEPhi_
    );

    phaseModel::update();
    thermo_.update();
}


void Foam::multiPhaseModel::calcAlphaAndRho()
{
    alphaRho_ = dimensionedScalar("0", dimDensity, 0.0);

    // find largest volume fraction and set to 1-sum
    forAll(rho_, celli)
    {
        SortableList<scalar> alphas(alphas_.size());
        forAll(alphas_, phasei)
        {
            alphas_[phasei][celli] =
                Foam::max
                (
                    Foam::min
                    (
                        alphas_[phasei][celli],
                        1.0
                    ),
                    0.0
                );
            alphas[phasei] = alphas_[phasei][celli];
        }
        alphas.reverseSort();

        const label fixedPhase = alphas.indices()[0];

        scalar sumAlpha = 0.0;
        for (label phasei = 1; phasei < alphas.size(); phasei++)
        {
            sumAlpha += alphas[phasei];
        }
        alphas_[fixedPhase][celli] = (*this)[celli] - sumAlpha;
    }

    forAll(alphas_, phasei)
    {
        alphas_[phasei].correctBoundaryConditions();

        alphaRhos_[phasei].max(0);
        rhos_[phasei] =
            alphaRhos_[phasei]
           /Foam::max
            (
                alphas_[phasei],
                thermo_.thermo(phasei).residualAlpha()
            );
        rhos_[phasei].correctBoundaryConditions();

        alphaRhos_[phasei] = alphas_[phasei]*rhos_[phasei];

        alphaRho_ += alphaRhos_[phasei];
    }
}

void Foam::multiPhaseModel::decode()
{
    calcAlphaAndRho();
    volScalarField& alpha(*this);
    this->correctBoundaryConditions();

    rho_ = alphaRho_/Foam::max(alpha, residualAlpha());

    volScalarField alphaRho(alphaRho_);
    alphaRho.max(1e-10);
    U_.ref() = alphaRhoU_()/(alphaRho());
    U_.correctBoundaryConditions();

    alphaRhoU_.boundaryFieldRef() =
        alphaRho_.boundaryField()*U_.boundaryField();

    volScalarField E(alphaRhoE_/alphaRho);
    e_.ref() = E() - 0.5*magSqr(U_());

    //- Limit internal energy it there is a negative temperature
    if (Foam::min(this->T_).value() < 0.0 && thermo_.limit())
    {
        if (debug)
        {
            WarningInFunction
                << "Lower limit of temperature reached, min(T) = "
                << Foam::min(T_).value()
                << ", limiting internal energy." << endl;
        }
        volScalarField limit(pos(T_));
        T_.max(small);
        e_ = e_*limit + thermo_.E()*(1.0 - limit);
        alphaRhoE_.ref() = rho_*(e_() + 0.5*magSqr(U_()));
    }
    e_.correctBoundaryConditions();

    alphaRhoE_.boundaryFieldRef() =
        alphaRho_.boundaryField()
       *(
            e_.boundaryField()
          + 0.5*magSqr(U_.boundaryField())
        );

    thermo_.correct();
}


void Foam::multiPhaseModel::encode()
{

    //- Scale volume fractions to new value
    alphaRho_ = dimensionedScalar(dimDensity, 0.0);
    volScalarField& alpha(*this);
    alpha = 0.0;
    forAll(alphas_, phasei)
    {
        alpha += alphas_[phasei];
        alphaRhos_[phasei] = alphas_[phasei]*rhos_[phasei];
        alphaRho_ += alphaRhos_[phasei];
    }

    alphaRhoU_ = alphaRho_*U_;
    alphaRhoE_ = alphaRho_*(e_ + 0.5*magSqr(U_));
}


Foam::tmp<Foam::volVectorField>
Foam::multiPhaseModel::gradP() const
{
    return fvc::grad(fluxScheme_->pf());
}


Foam::tmp<Foam::volVectorField>
Foam::multiPhaseModel::gradAlpha() const
{
    return fvc::grad(fluxScheme_->alphaf());
}



void Foam::multiPhaseModel::correctThermo()
{
    thermo_.correct();
}


Foam::tmp<Foam::volScalarField>
Foam::multiPhaseModel::speedOfSound() const
{
    return thermo_.speedOfSound();
}


Foam::tmp<Foam::volScalarField>
Foam::multiPhaseModel::ESource() const
{
    return thermo_.ESource();
}


Foam::tmp<Foam::volScalarField> Foam::multiPhaseModel::Cv() const
{
    return thermo_.Cv()/Foam::max(*this, residualAlpha());
}

Foam::tmp<Foam::volScalarField> Foam::multiPhaseModel::Cp() const
{
    return thermo_.Cp()/Foam::max(*this, residualAlpha());
}

Foam::scalar Foam::multiPhaseModel::Cvi(const label celli) const
{
    return
        thermo_.Cvi(celli)
       /Foam::max(this->operator[](celli), residualAlpha().value());
}

Foam::tmp<Foam::volScalarField> Foam::multiPhaseModel::mu() const
{
    return thermo_.mu()/Foam::max(*this, residualAlpha());
}


Foam::tmp<Foam::scalarField>
Foam::multiPhaseModel::mu(const label patchi) const
{
    return
        thermo_.mu(patchi)
       /Foam::max(this->boundaryField()[patchi], residualAlpha().value());
}

Foam::tmp<Foam::volScalarField> Foam::multiPhaseModel::nu() const
{
    return thermo_.nu()/Foam::max(*this, residualAlpha());
}

Foam::tmp<Foam::scalarField>
Foam::multiPhaseModel::nu(const label patchi) const
{
    return
        thermo_.nu(patchi)
       /Foam::max(this->boundaryField()[patchi], residualAlpha().value());
}

Foam::scalar
Foam::multiPhaseModel::nui(const label celli) const
{
    return
        thermo_.nui(celli)
       /Foam::max(this->operator[](celli), residualAlpha().value());
}

Foam::tmp<Foam::volScalarField> Foam::multiPhaseModel::alpha() const
{
    return thermo_.alpha()/Foam::max(*this, residualAlpha());
}

Foam::tmp<Foam::scalarField>
Foam::multiPhaseModel::alpha(const label patchi) const
{
    return
        thermo_.alpha(patchi)
       /Foam::max(this->boundaryField()[patchi], residualAlpha().value());
}

Foam::tmp<Foam::volScalarField> Foam::multiPhaseModel::alphaEff
(
    const volScalarField& alphat
) const
{
    return thermo_.alphaEff(alphat)/Foam::max(*this, residualAlpha());
}

Foam::tmp<Foam::scalarField> Foam::multiPhaseModel::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        thermo_.alphaEff(alphat, patchi)
       /Foam::max(this->boundaryField()[patchi], residualAlpha().value());
}

Foam::tmp<Foam::volScalarField> Foam::multiPhaseModel::alphahe() const
{
    return thermo_.alphahe()/Foam::max(*this, residualAlpha());
}

Foam::tmp<Foam::scalarField>
Foam::multiPhaseModel::alphahe(const label patchi) const
{
    return
        thermo_.alphahe(patchi)
       /Foam::max(this->boundaryField()[patchi], residualAlpha().value());
}

Foam::tmp<Foam::volScalarField> Foam::multiPhaseModel::kappa() const
{
    return thermo_.kappa()/Foam::max(*this, residualAlpha());
}

Foam::tmp<Foam::scalarField>
Foam::multiPhaseModel::kappa(const label patchi) const
{
    return
        thermo_.kappa(patchi)
       /Foam::max(this->boundaryField()[patchi], residualAlpha().value());
}

Foam::scalar
Foam::multiPhaseModel::kappai(const label celli) const
{
    return
        thermo_.kappai(celli)
       /Foam::max(this->operator[](celli), residualAlpha().value());
}

Foam::tmp<Foam::volScalarField> Foam::multiPhaseModel::kappaEff
(
    const volScalarField& alphat
) const
{
    return thermo_.kappaEff(alphat)/Foam::max(*this, residualAlpha());
}

Foam::tmp<Foam::scalarField>
Foam::multiPhaseModel::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        thermo_.kappaEff(alphat, patchi)
       /Foam::max(this->boundaryField()[patchi], residualAlpha().value());
}

// ************************************************************************* //
