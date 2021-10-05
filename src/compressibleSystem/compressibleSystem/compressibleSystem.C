/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "compressibleSystem.H"
#include "uniformDimensionedFields.H"
#include "fvm.H"
#include "wedgeFvPatch.H"
#include "blastRadiationModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(compressibleSystem, 0);
    defineRunTimeSelectionTable(compressibleSystem, dictionary);
}


// * * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * //

void Foam::compressibleSystem::setModels()
{
    if (Foam::max(this->thermo().mu()).value() > small)
    {
        needPostUpdate_ = true;
        turbulence_ =
        (
            compressible::momentumTransportModel::New
            (
                rho(),
                U(),
                rhoPhi(),
                this->thermo()
            )
        );
        turbulence_->validate();

        thermophysicalTransport_ =
        (
            fluidThermophysicalTransportModel::New
            (
                turbulence_,
                this->thermo()
            ).ptr()
        );
    }


    modelsPtr_.set(&fvModels::New(mesh()));
    constraintsPtr_.set(&fvConstraints::New(mesh()));

    const PtrList<fvModel>& models(modelsPtr_());
    forAll(models, modeli)
    {
        wordList fields(models[modeli].addSupFields());
        forAll(fields, fieldi)
        {
            if (!solveFields_.found(fields[fieldi]))
            {
                solveFields_.append(fields[fieldi]);
            }
        }
    }

    const PtrList<fvConstraint>& constraints(constraintsPtr_());
    forAll(constraints, modeli)
    {
        wordList fields(constraints[modeli].constrainedFields());
        forAll(fields, fieldi)
        {
            if (!solveFields_.found(fields[fieldi]))
            {
                solveFields_.append(fields[fieldi]);
            }
        }
    }
    needPostUpdate_ = solveFields_.size();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressibleSystem::compressibleSystem
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "phaseProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    timeIntegrationSystem("compressibleSystem", mesh),
    U_
    (
        IOobject
        (
            "U",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    rhoU_
    (
        IOobject
        (
            "rhoU",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimDensity*dimVelocity, Zero),
        "zeroGradient"
    ),
    rhoE_
    (
        IOobject
        (
            "rhoE",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimDensity*sqr(dimVelocity), 0.0)
    ),
    phi_
    (
        IOobject
        (
            "phi",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimVelocity*dimArea, 0.0)
    ),
    rhoPhi_
    (
        IOobject
        (
            "rhoPhi",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimDensity*dimVelocity*dimArea, 0.0)
    ),
    rhoUPhi_
    (
        IOobject
        (
            "rhoUPhi",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedVector("0", dimDensity*sqr(dimVelocity)*dimArea, Zero)
    ),
    rhoEPhi_
    (
        IOobject
        (
            "rhoEPhi",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimDensity*pow3(dimVelocity)*dimArea, 0.0)
    ),
    fluxScheme_(fluxScheme::New(mesh)),
    g_(mesh.lookupObject<uniformDimensionedVectorField>("g")),
    solutionDs_((vector(mesh.solutionD()) + vector::one)/2.0),
    solveFields_(),
    needPostUpdate_(false)
{
    scalar emptyDirV
    (
        Foam::max(mag(U_ & (vector::one - solutionDs_))).value()
    );

    // Remove wedge directions if not used
    if (emptyDirV < small)
    {
        solutionDs_ = ((vector(mesh.geometricD()) + vector::one)/2.0);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::compressibleSystem::~compressibleSystem()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::compressibleSystem::encode()
{
    rhoU_ = rho()*U_;
    rhoE_ = rho()*(he() + 0.5*magSqr(U_));
}


void Foam::compressibleSystem::update()
{
    decode();
    fluxScheme_->update
    (
        rho(),
        U(),
        he(),
        p(),
        speedOfSound()(),
        phi_,
        rhoPhi_,
        rhoUPhi_,
        rhoEPhi_
    );
}


void Foam::compressibleSystem::solve()
{
    //- Calculate deltas for momentum and energy
    volVectorField deltaRhoU
    (
        "deltaRhoU",
        fvc::div(rhoUPhi_) - g_*rho()
    );

    volScalarField deltaRhoE
    (
        "deltaRhoE",
        fvc::div(rhoEPhi_)
      - (rhoU_ & g_)
    );

    //- Store old values
    this->storeAndBlendOld(rhoU_);
    this->storeAndBlendOld(rhoE_);

    //- Store changed in momentum and energy
    this->storeAndBlendDelta(deltaRhoU);
    this->storeAndBlendDelta(deltaRhoE);

    //- Solve for momentum and energy
    dimensionedScalar dT = rho().time().deltaT();
    rhoU_ -= cmptMultiply(dT*deltaRhoU, solutionDs_);
    rhoE_ -= dT*deltaRhoE;
}


void Foam::compressibleSystem::postUpdate()
{
    // Solve momentum
    if (solveFields_.found(U_.name()) || turbulence_.valid())
    {
        fvVectorMatrix UEqn
        (
            fvm::ddt(rho(), U_) - fvc::ddt(rho(), U_)
        ==
            modelsPtr_->source(rho(), U_)
        );
        if (turbulence_.valid())
        {
            UEqn += turbulence_->divDevTau(U_);
            rhoE_ +=
                rho().mesh().time().deltaT()
                *fvc::div
                (
                    fvc::dotInterpolate(rho().mesh().Sf(), turbulence_->devTau())
                    & fluxScheme_->Uf()
                );
        }
        constraintsPtr_->constrain(UEqn);
        UEqn.solve();
        constraintsPtr_->constrain(U_);

        //- Update internal energy
        he() = rhoE_/rho() - 0.5*magSqr(U_);
    }

    // Solve thermal energy diffusion
    if (solveFields_.found(he().name()) || turbulence_.valid())
    {
        fvScalarMatrix eEqn
        (
            fvm::ddt(rho(), he()) - fvc::ddt(rho(), he())
        ==
            modelsPtr_->source(rho(), he())
        );
        if (turbulence_.valid())
        {
            eEqn += thermophysicalTransport_->divq(he());
        }
        constraintsPtr_->constrain(eEqn);
        eEqn.solve();
        constraintsPtr_->constrain(he());
    }

    if (turbulence_.valid())
    {
        turbulence_->correct();
    }

    encode();
    this->thermo().correct();
}


Foam::scalar Foam::compressibleSystem::CoNum() const
{
    surfaceScalarField amaxSf
    (
        fvc::interpolate(speedOfSound())*mesh().magSf()
    );
    // Remove wave speed from wedge boundaries
    forAll(amaxSf.boundaryField(), patchi)
    {
        if (isA<wedgeFvPatch>(mesh().boundary()[patchi]))
        {
            amaxSf.boundaryFieldRef() = Zero;
        }
    }
    amaxSf += mag(fvc::flux(U()));

    scalarField sumAmaxSf
    (
        fvc::surfaceSum(amaxSf)().primitiveField()
    );

    return 0.5*gMax(sumAmaxSf/mesh().V().field())*mesh().time().deltaTValue();
}


bool Foam::compressibleSystem::writeData(Ostream& os) const
{
    return os.good();
}


bool Foam::compressibleSystem::read()
{
    return regIOobject::read();
}

// ************************************************************************* //
