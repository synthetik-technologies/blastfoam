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
#include "phaseFluxScheme.H"
#include "interfacialPressureModel.H"
#include "interfacialVelocityModel.H"

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
    timeIntegrationSystem
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
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh()
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
    solutionDs_((vector(fluid.mesh().solutionD()) + vector::one)/2.0)

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
    dPtr_.reset
    (
        diameterModel::New(fluid_.mesh(), phaseDict_, name_).ptr()
    );
    thermo().initializeModels();
    thermo().correct();
}


Foam::autoPtr<Foam::phaseModel> Foam::phaseModel::clone() const
{
    NotImplemented;
    return autoPtr<phaseModel>(nullptr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseModel::~phaseModel()
{
    if (dPtr_.valid())
    {
        delete dPtr_.ptr();
    }
}


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
    volScalarField PI
    (
        volScalarField::New
        (
            IOobject::groupName("PI", name_),
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
            PI +=
                fluid_.lookupSubModel<interfacialPressureModel>
                (
                    *this,
                    otherPhase,
                    false
                ).PI()
               *otherPhase;
            sumAlpha.ref() += otherPhase;
        }
    }

    PI /= Foam::max(sumAlpha, this->residualAlpha());
    dPtr_->solve(PI, T());
}


void Foam::phaseModel::solveAlphaRho()
{
    volScalarField deltaAlphaRho(fvc::div(alphaRhoPhi_));
    forAll(fluid_.phases(), phasei)
    {
        const phaseModel& otherPhase = fluid_.phases()[phasei];
        if (&otherPhase != this)
        {
            deltaAlphaRho -= fluid_.mDot(*this, otherPhase);
        }
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

    tmp<volScalarField> deltaAlpha;
    if (solveAlpha_)
    {
        volScalarField& alpha(*this);
        deltaAlpha =
            volScalarField::New
            (
                IOobject::groupName("deltaAlpha", name_),
                fvc::div(alphaPhiPtr_()) - alpha*fvc::div(fluid_.phi())
            );
    }

    volVectorField deltaAlphaRhoU
    (
        IOobject::groupName("deltaAlphaRhoU", name_),
        fvc::div(alphaRhoUPhi_)
      - fluid_.PI()*gradAlpha()
      - (*this)*rho()*fluid_.g() // alphaRho has already been updated
    );

    volScalarField deltaAlphaRhoE
    (
        IOobject::groupName("deltaAlphaRhoE", name_),
        fvc::div(alphaRhoEPhi_)
      - ESource()
      - fluid_.PI()*(fluid_.U() & gradAlpha())
      - (alphaRhoU_ & fluid_.g())
    );

    forAll(fluid_.phases(), phasei)
    {
        const phaseModel& otherPhase = fluid_.phases()[phasei];
        if (&otherPhase != this)
        {
            volScalarField mD(fluid_.mDot(*this, otherPhase));
            if (solveAlpha_)
            {
                deltaAlpha.ref() -=
                    fluid_.mDotByRho(mD, *this, otherPhase);
            }
            deltaAlphaRhoU -= fluid_.mDotU(mD, *this, otherPhase);
            deltaAlphaRhoE -= fluid_.mDotE(mD, *this, otherPhase);
        }
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
        this->storeAndBlendDelta(deltaAlpha.ref());

        volScalarField& alpha(*this);
        this->storeAndBlendOld(alpha);

        alpha -= dT*deltaAlpha;
        alpha.max(0);
        alpha.min(alphaMax_);
        alpha.correctBoundaryConditions();
    }
}


void Foam::phaseModel::postUpdate()
{
    dimensionedScalar smallAlphaRho(residualAlphaRho());
    if (needSolve(U_.name()) || turbulence_.valid())
    {
        fvVectorMatrix UEqn
        (
            fvm::ddt(alphaRho_, U_) - fvc::ddt(alphaRhoU_)
          + fvc::ddt(smallAlphaRho, U_) - fvm::ddt(smallAlphaRho, U_)
         ==
            models().source(alphaRho_, U_)
        );
        if (turbulence_.valid())
        {
            UEqn += turbulence_->divDevTau(U_);
            alphaRhoE_ +=
                rho().time().deltaT()
               *fvc::div
                (
                    fvc::dotInterpolate
                    (
                        rho().mesh().Sf(),
                        turbulence_->devTau()
                    )
                  & flux().Uf()
                );
        }
        constraints().constrain(UEqn);
        UEqn.solve();
        constraints().constrain(U_);

        alphaRhoU_ = alphaRho_*U_;

        he() = alphaRhoE_/Foam::max(alphaRho_, smallAlphaRho) - 0.5*magSqr(U_);
    }

    // Solve thermal energy diffusion
    if (needSolve(he().name()) || turbulence_.valid())
    {
        fvScalarMatrix eEqn
        (
            fvm::ddt(alphaRho_, he())
          - fvc::ddt(alphaRho_.prevIter(), he())
//           - fvc::ddt(alphaRhoE_)
//           - (
//                 U_
//               & (
//                     alphaRho_*fvc::ddt(U_)
//                   + 0.5*U_*fvc::ddt(alphaRho_)
//                 )
//             )
          + fvc::ddt(smallAlphaRho, he())
          - fvm::ddt(smallAlphaRho, he())
         ==
            models().source(alphaRho_, he())
        );

        if (turbulence_.valid())
        {
            // Add thermal energy diffusion
            eEqn += thermophysicalTransport_->divq(he());
        }
        constraints().constrain(eEqn);
        eEqn.solve();
        constraints().constrain(he());

        alphaRhoE_ = alphaRho_*(he() + 0.5*magSqr(U_));
    }

    if (turbulence_.valid())
    {
        turbulence_->correct();
    }

    thermo().postUpdate();
    thermo().correct();
}


void Foam::phaseModel::correctVolumeFraction()
{}


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
