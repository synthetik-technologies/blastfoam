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
            fluid.mesh().time().name(),
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
            fluid.mesh().time().name(),
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
            fluid.mesh().time().name(),
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
            fluid.mesh().time().name(),
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
            fluid.mesh().time().name(),
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
            fluid.mesh().time().name(),
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
            fluid.mesh().time().name(),
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
            fluid.mesh().time().name(),
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
}


void Foam::phaseModel::solveD()
{
    if (!dPtr_->requireInput())
    {
        dPtr_->solve();
        return;
    }

    // Evolve the diameter model using the interfacial pressure and the
    // surface temperature
    dPtr_->solve(fluid_.p(), Ts());
}


void Foam::phaseModel::solveAlphaRho()
{
    volScalarField deltaAlphaRho
    (
        IOobject::groupName("deltaAlphaRho", name_),
        fvc::div(alphaRhoPhi_)
    );
    this->fvTimeInt_->addDeltaSource(alphaRho_.name(), deltaAlphaRho);
    if (fluid_.hasMassTransfer(*this))
    {
        forAll(fluid_.phases(), phasei)
        {
            const phaseModel& otherPhase = fluid_.phases()[phasei];
            if (&otherPhase != this)
            {
                deltaAlphaRho -= fluid_.mDot(*this, otherPhase);
            }
        }
    }

    this->storeAndBlendDelta(deltaAlphaRho);


    this->storeAndBlendOld(alphaRho_);
    alphaRho_.storePrevIter();
    alphaRho_ -= this->mesh().time().deltaT()*deltaAlphaRho;
    alphaRho_.max(0);
}


void Foam::phaseModel::scaleVolumeFraction
(
    const scalar sumAlpha,
    const label celli
)
{
    (*this)[celli] /= sumAlpha;
}


void Foam::phaseModel::correctVolumeFraction
(
    const scalar alpha,
    const label celli
)
{
    (*this)[celli] = alpha;
}


void Foam::phaseModel::update()
{
    dPtr_->update();
    solveD();
}


Foam::tmp<Foam::surfaceScalarField> Foam::phaseModel::alphaPhi() const
{
    return this->fluid_.phi()*flux().alphaf();
}


// ************************************************************************* //
