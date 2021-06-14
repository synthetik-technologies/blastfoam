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

#include "phaseModel.H"
#include "phaseSystem.H"
#include "fvMatrix.H"
#include "fixedValueFvsPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"
#include "basicThermoModel.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "fluxScheme.H"
#include "interfacialPressureModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseModel, 0);
    defineRunTimeSelectionTable(phaseModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseModel::phaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const label index
)
:
    volScalarField
    (
        IOobject
        (
            IOobject::groupName("alpha", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("alpha", dimless, 0)
    ),
    integrationSystem
    (
        IOobject::groupName("phaseModel", phaseName),
        fluid.mesh()
    ),
    fluid_(fluid),
    name_(phaseName),
    index_(index),
    phaseDict_
    (
        fluid_.subDict(name_)
    ),
    alphaMax_(phaseDict_.lookupOrDefault("alphaMax", 1.0)),
    U_
    (
        IOobject
        (
            IOobject::groupName("U", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedVector("0", dimVelocity, Zero)
    ),
    T_
    (
        IOobject
        (
            IOobject::groupName("T", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh()
    ),
    e_
    (
        IOobject
        (
            IOobject::groupName("e", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fluid_.mesh(),
        dimensionedScalar(sqr(dimVelocity), -1.0),
        basicThermoModel::eBoundaryTypes(T_),
        basicThermoModel::eBoundaryBaseTypes(T_)
    ),
    alphaRho_
    (
        IOobject
        (
            IOobject::groupName("alphaRho", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar("0", dimDensity, 0.0)
    ),
    alphaRhoU_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoU", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedVector("0", dimDensity*dimVelocity, Zero),
        "zeroGradient"
    ),
    alphaRhoE_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoE", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar("0", dimDensity*sqr(dimVelocity), 0.0)
    ),
    phi_(IOobject::groupName("phi", name_), fvc::flux(U_)),
    alphaRhoPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoPhi", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar("0", dimDensity*phi_.dimensions(), 0.0)
    ),
    alphaRhoUPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoUPhi", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedVector("0", dimVelocity*alphaRhoPhi_.dimensions(), Zero)
    ),
    alphaRhoEPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoEPhi", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar("0", dimVelocity*alphaRhoUPhi_.dimensions(), 0.0)
    ),
    solutionDs_((vector(this->mesh().solutionD()) + vector::one)/2.0)
{
    scalar emptyDirV
    (
        Foam::max(mag(U_ & (vector::one - solutionDs_))).value()
    );

    // Remove wedge directions if not used
    if (emptyDirV < small)
    {
        solutionDs_ = ((vector(this->mesh().geometricD()) + vector::one)/2.0);
    }
}


void Foam::phaseModel::initializeModels()
{
    dPtr_ = diameterModel::New(fluid_.mesh(), phaseDict_, name_);
}


Foam::autoPtr<Foam::phaseModel> Foam::phaseModel::clone() const
{
    NotImplemented;
    return autoPtr<phaseModel>(nullptr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseModel::~phaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::phaseModel::solveAlpha(const bool s)
{
    solveAlpha_ = s;
    if (alphaPhiPtr_.valid() || !s)
    {
        return;
    }

    alphaPhiPtr_.set
    (
        new surfaceScalarField
        (
            IOobject
            (
                IOobject::groupName("alphaPhi", name_),
                this->mesh().time().timeName(),
                this->mesh()
            ),
            this->mesh(),
            dimensionedScalar("0", phi_.dimensions(), 0.0)
        )
    );
}


void Foam::phaseModel::solveD()
{
    if (!dPtr_->requireInput())
    {
        dPtr_->solve();
        return;
    }

    //- Update diameterModel
    volScalarField Pi
    (
        volScalarField::New
        (
            IOobject::groupName("pi", name_),
            fluid_.mesh(),
            dimensionedScalar(dimPressure, 0.0)
        )
    );
    tmp<volScalarField> sumAlpha
    (
        volScalarField::New
        (
            IOobject::groupName("sumAlpha", name_),
            fluid_.mesh(),
            dimensionedScalar(dimless, 0.0)
        )
    );

    forAll(fluid_.phases(), phasei)
    {
        const phaseModel& otherPhase = fluid_.phases()[phasei];
        if
        (
            fluid_.foundSubModel<interfacialPressureModel>
            (
                *this,
                otherPhase,
                false
            )
        )
        {
            Pi +=
                fluid_.lookupSubModel<interfacialPressureModel>
                (
                    *this,
                    otherPhase,
                    false
                ).Pi()
               *otherPhase;
            sumAlpha.ref() += otherPhase;
        }
    }

    Pi /= Foam::max(sumAlpha, this->residualAlpha());
    dPtr_->solve(Pi, T_);
}


void Foam::phaseModel::solveAlphaRho()
{
    volScalarField deltaAlphaRho(fvc::div(alphaRhoPhi_));
    forAll(fluid_.phases(), phasei)
    {
        const phaseModel& phase = fluid_.phases()[phasei];
        deltaAlphaRho -= fluid_.mDot(*this, phase);
    }

    this->storeAndBlendDelta(deltaAlphaRho);

    this->storeAndBlendOld(alphaRho_);
    alphaRho_.storePrevIter();
    alphaRho_ -= this->mesh().time().deltaT()*deltaAlphaRho;
    alphaRho_.max(0);
    alphaRho_.correctBoundaryConditions();
}


void Foam::phaseModel::solve()
{
    dimensionedScalar dT = rho().time().deltaT();

    volVectorField deltaAlphaRhoU
    (
        IOobject::groupName("deltaAlphaRhoU", name_),
        fvc::div(alphaRhoUPhi_)
      - p()*gradAlpha()
      - (*this)*rho()*fluid_.g() // alphaRho has already been updated
    );

    volScalarField deltaAlphaRhoE
    (
        IOobject::groupName("deltaAlphaRhoE", name_),
        fvc::div(alphaRhoEPhi_)
      - ESource()
      - (alphaRhoU_ & fluid_.g())
    );

    forAll(fluid_.phases(), phasei)
    {
        const phaseModel& phase = fluid_.phases()[phasei];
        if (phase.granular())
        {
            deltaAlphaRhoE += p()*fvc::div(phase.alphaPhi());
        }
        deltaAlphaRhoU -= fluid_.mDotU(*this, phase);
        deltaAlphaRhoE -= fluid_.mDotE(*this, phase);
    }
    this->storeAndBlendDelta(deltaAlphaRhoU);
    this->storeAndBlendDelta(deltaAlphaRhoE);


    this->storeAndBlendOld(alphaRhoU_);
    alphaRhoU_ -= cmptMultiply(dT*deltaAlphaRhoU, solutionDs_);
    alphaRhoU_.correctBoundaryConditions();

    this->storeAndBlendOld(alphaRhoE_);
    alphaRhoE_ -= dT*deltaAlphaRhoE;
    alphaRhoE_.correctBoundaryConditions();

    // Transport volume fraction if required
    if (solveAlpha_)
    {
        volScalarField& alpha = *this;
        volScalarField deltaAlpha
        (
            IOobject::groupName("deltaAlpha", name_),
            fvc::div(alphaPhiPtr_()) - alpha*fvc::div(fluid_.phi())
        );
        this->storeAndBlendDelta(deltaAlpha);

        this->storeAndBlendOld(alpha);
        alpha -= dT*deltaAlpha;
        alpha.max(0);
        alpha.min(alphaMax_);
        alpha.correctBoundaryConditions();
    }
}


void Foam::phaseModel::postUpdate()
{
    if (turbulence_.valid())
    {
        dimensionedScalar smallAlphaRho(dimDensity, 1e-6);
        fvVectorMatrix UEqn
        (
            fvm::ddt(alphaRho_, U_) - fvc::ddt(alphaRho_, U_)
          + fvc::ddt(smallAlphaRho, U_) - fvm::ddt(smallAlphaRho, U_)
          + turbulence_->divDevTau(U_)
        );

        alphaRhoE_ +=
            rho().time().deltaT()
           *fvc::div
            (
                fvc::dotInterpolate(rho().mesh().Sf(), turbulence_->devTau())
              & flux().Uf()
            );

        UEqn.solve();
        alphaRhoU_ = cmptMultiply(alphaRho_*U_, solutionDs_);

        // Solve thermal energy diffusion
        e_ = alphaRhoE_/Foam::max(alphaRho_, smallAlphaRho) - 0.5*magSqr(U_);
        fvScalarMatrix eEqn
        (
            fvm::ddt(alphaRho_, e_) - fvc::ddt(alphaRho_, e_)
          + fvc::ddt(smallAlphaRho, e_) - fvm::ddt(smallAlphaRho, e_)
          + thermophysicalTransport_->divq(e_)
        );
        eEqn.solve();
        alphaRhoE_ = alphaRho_*(e_ + 0.5*magSqr(U_));

        turbulence_->correct();
    }
}


void Foam::phaseModel::update()
{
    dPtr_->update();
    solveD();
}


Foam::tmp<Foam::surfaceScalarField> Foam::phaseModel::alphaPhi() const
{
    if (alphaPhiPtr_.valid())
    {
        return alphaPhiPtr_();
    }
    return this->phi_*fvc::interpolate(*this);
}


// ************************************************************************* //
