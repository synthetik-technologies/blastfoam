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
#include "PhaseCompressibleTurbulenceModel.H"

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
        fluid.mesh(),
        dimensionedScalar("0", sqr(dimVelocity), -1),
        T_.boundaryField().types()
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
        dimensionedVector("0", dimDensity*dimVelocity, Zero)
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
    dPtr_(diameterModel::New(fluid.mesh(), phaseDict_, phaseName))
{}


Foam::autoPtr<Foam::phaseModel> Foam::phaseModel::clone() const
{
    NotImplemented;
    return autoPtr<phaseModel>(nullptr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseModel::~phaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::phaseModel::solve
(
    const label stepi,
    const scalarList& ai,
    const scalarList& bi
)
{
    if (oldIs_[stepi - 1] != -1)
    {
        alphaRhoOld_.set(oldIs_[stepi - 1], new volScalarField(alphaRho_));
        alphaRhoUOld_.set
        (
            oldIs_[stepi - 1],
            new volVectorField(alphaRhoU_)
        );
        alphaRhoEOld_.set
        (
            oldIs_[stepi - 1],
            new volScalarField(alphaRhoE_)
        );
    }
    volScalarField alphaRhoOld(ai[stepi - 1]*alphaRho_);
    volVectorField alphaRhoUOld(ai[stepi - 1]*alphaRhoU_);
    volScalarField alphaRhoEOld(ai[stepi - 1]*alphaRhoE_);

    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = oldIs_[i];
        if (fi != -1 && ai[fi] != 0)
        {
            alphaRhoOld += ai[fi]*alphaRhoOld_[fi];
            alphaRhoUOld += ai[fi]*alphaRhoUOld_[fi];
            alphaRhoEOld += ai[fi]*alphaRhoEOld_[fi];
        }
    }

    volScalarField deltaAlphaRho(fvc::div(alphaRhoPhi_));
    volVectorField deltaAlphaRhoU
    (
        fvc::div(alphaRhoUPhi_)
      - p()*gradAlpha()
    );
    volScalarField deltaAlphaRhoE
    (
        fvc::div(alphaRhoEPhi_) - ESource()
    );
    forAll(fluid_.phases(), phasei)
    {
        const phaseModel& phase = fluid_.phases()[phasei];
        if (phase.granular())
        {
            deltaAlphaRhoE += p()*fvc::div(phase.alphaPhi());
        }
    }

    if (deltaIs_[stepi - 1] != -1)
    {
        deltaAlphaRho_.set
        (
            deltaIs_[stepi - 1],
            new volScalarField(deltaAlphaRho)
        );
        deltaAlphaRhoU_.set
        (
            deltaIs_[stepi - 1],
            new volVectorField(deltaAlphaRhoU)
        );
        deltaAlphaRhoE_.set
        (
            deltaIs_[stepi - 1],
            new volScalarField(deltaAlphaRhoE)
        );
    }
    deltaAlphaRho *= bi[stepi - 1];
    deltaAlphaRhoU *= bi[stepi - 1];
    deltaAlphaRhoE *= bi[stepi - 1];

    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = deltaIs_[i];
        if (fi != -1 && bi[fi] != 0)
        {
            deltaAlphaRho += bi[fi]*deltaAlphaRho_[fi];
            deltaAlphaRhoU += bi[fi]*deltaAlphaRhoU_[fi];
            deltaAlphaRhoE += bi[fi]*deltaAlphaRhoE_[fi];
        }
    }


    dimensionedScalar dT = this->mesh().time().deltaT();
    vector solutionDs((vector(this->mesh().solutionD()) + vector::one)/2.0);

    alphaRho_ = alphaRhoOld - dT*(deltaAlphaRho);
    alphaRho_.max(0);
    alphaRhoU_ = cmptMultiply(alphaRhoUOld - dT*deltaAlphaRhoU, solutionDs);
    alphaRhoE_ = alphaRhoEOld - dT*deltaAlphaRhoE;

    //- Update diameterModel
    dPtr_->solve(stepi, ai, bi);
}


void Foam::phaseModel::setODEFields
(
    const label nSteps,
    const boolList& storeFields,
    const boolList& storeDeltas
)
{
    oldIs_.resize(nSteps);
    deltaIs_.resize(nSteps);
    label fi = 0;
    for (label i = 0; i < nSteps; i++)
    {
        if (storeFields[i])
        {
            oldIs_[i] = fi;
            fi++;
        }
        else
        {
            oldIs_[i] = -1;
        }
    }
    nOld_ = fi;
    alphaRhoOld_.resize(nOld_);
    alphaRhoUOld_.resize(nOld_);
    alphaRhoEOld_.resize(nOld_);

    fi = 0;
    for (label i = 0; i < nSteps; i++)
    {
        if (storeDeltas[i])
        {
            deltaIs_[i] = fi;
            fi++;
        }
        else
        {
            deltaIs_[i] = -1;
        }
    }
    nDelta_ = fi;
    deltaAlphaRho_.resize(nDelta_);
    deltaAlphaRhoU_.resize(nDelta_);
    deltaAlphaRhoE_.resize(nDelta_);
}


void Foam::phaseModel::clearODEFields()
{
    alphaRhoOld_.clear();
    alphaRhoUOld_.clear();
    alphaRhoEOld_.clear();

    alphaRhoOld_.resize(nOld_);
    alphaRhoUOld_.resize(nOld_);
    alphaRhoEOld_.resize(nOld_);

    deltaAlphaRho_.clear();
    deltaAlphaRhoU_.clear();
    deltaAlphaRhoE_.clear();

    deltaAlphaRho_.resize(nDelta_);
    deltaAlphaRhoU_.resize(nDelta_);
    deltaAlphaRhoE_.resize(nDelta_);
}

// ************************************************************************* //
