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

#include "phaseCompressibleSystem.H"
#include "uniformDimensionedFields.H"
#include "fvm.H"
#include "wedgeFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseCompressibleSystem, 0);
    defineRunTimeSelectionTable(phaseCompressibleSystem, dictionary);
}


void Foam::phaseCompressibleSystem::setModels()
{
    if (Foam::max(this->thermo().mu()).value() > 0)
    {
        turbulence_ =
        (
            blast::momentumTransportModel::New
            (
                rho_,
                U_,
                rhoPhi_,
                this->thermo()
            )
        );
        turbulence_->validate();

        thermophysicalTransport_ =
        (
            blastFluidThermophysicalTransportModel::New
            (
                turbulence_,
                this->thermo()
            ).ptr()
        );
    }


    IOobject radPropertiesIO
    (
        "radiationProperties",
        rho_.time().constant(),
        rho_.mesh(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    );


    if (radPropertiesIO.typeHeaderOk<IOdictionary>(true))
    {
        radiation_.set(radiationModel::New(this->T()).ptr());
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseCompressibleSystem::phaseCompressibleSystem(const fvMesh& mesh)
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
    integrationSystem("phaseCompressibleSystem", mesh),
    rho_
    (
        IOobject
        (
            "rho",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rho", dimDensity, 0.0)
    ),
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
    p_
    (
        IOobject
        (
            "p",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    T_
    (
        IOobject
        (
            "T",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    e_
    (
        IOobject
        (
            "e",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(sqr(dimVelocity), -1.0),
        fluidThermoModel::eBoundaryTypes(T_),
        fluidThermoModel::eBoundaryBaseTypes(T_)
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
        rho_*U_,
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
    TLow_("TLow", dimTemperature, 0.0),
    solutionDs_((vector(mesh.solutionD()) + vector::one)/2.0)
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

    TLow_.readIfPresent(*this);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseCompressibleSystem::~phaseCompressibleSystem()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::phaseCompressibleSystem::solve()
{
    //- Calculate deltas for momentum and energy
    volVectorField deltaRhoU
    (
        "deltaRhoU",
        fvc::div(rhoUPhi_) - g_*rho_
    );

    volScalarField deltaRhoE
    (
        "deltaRhoE",
        fvc::div(rhoEPhi_)
      - ESource()
      - (rhoU_ & g_)
    );

    //- Store old values
    this->storeAndBlendOld(rhoU_);
    this->storeAndBlendOld(rhoE_);

    //- Store changed in momentum and energy
    this->storeAndBlendDelta(deltaRhoU);
    this->storeAndBlendDelta(deltaRhoE);

    //- Solve for momentum and energy
    dimensionedScalar dT = rho_.time().deltaT();
    rhoU_ -= cmptMultiply(dT*deltaRhoU, solutionDs_);
    rhoE_ -= dT*deltaRhoE;
}


void Foam::phaseCompressibleSystem::postUpdate()
{
    this->decode();

    this->thermo().postUpdate();

    if (radiation_.valid())
    {
        radiation_->correct();
        rhoE_ =
            radiation_->calcRhoE
            (
                rho_.mesh().time().deltaT(),
                rhoE_,
                rho_,
                e_,
                this->thermo().Cv()
            );
    }

    if
    (
        dragSource_.valid()
     || extESource_.valid()
     || turbulence_.valid()
    )
    {
        // Solve momentum
        fvVectorMatrix UEqn
        (
            fvm::ddt(rho_, U_) - fvc::ddt(rho_, U_)
        );

        if (dragSource_.valid())
        {
            UEqn -= dragSource_;
        }

        if (turbulence_.valid())
        {
            UEqn += turbulence_->divDevTau(U_);
            rhoE_ +=
                rho_.mesh().time().deltaT()
               *fvc::div
                (
                    fvc::dotInterpolate(rho_.mesh().Sf(), turbulence_->devTau())
                  & fluxScheme_->Uf()
                );
        }
        UEqn.solve();
        rhoU_ = cmptMultiply(rho_*U_, solutionDs_);

        // Solve thermal energy diffusion
        e_ = rhoE_/rho_ - 0.5*magSqr(U_);
        fvScalarMatrix eEqn
        (
            fvm::ddt(rho_, e_) - fvc::ddt(rho_, e_)
        );
        if (extESource_.valid())
        {
            eEqn -= extESource_;
        }
        if (turbulence_.valid())
        {
            eEqn += thermophysicalTransport_->divq(e_);
        }
        eEqn.solve();

        rhoE_ = rho_*(e_ + 0.5*magSqr(U_));
    }

    if (turbulence_.valid())
    {
        turbulence_->correct();
    }

    this->thermo().correct();
}


void Foam::phaseCompressibleSystem::addECoeff
(
    const volScalarField::Internal& ECoeff
)
{
    if (!extESource_.valid())
    {
        extESource_ =
            tmp<fvScalarMatrix>
            (
                new fvScalarMatrix(e(), dimEnergy/dimTime)
            );
    }
    extESource_.ref() -= fvm::Sp(ECoeff, e());
}


void Foam::phaseCompressibleSystem::addESource
(
    const volScalarField::Internal& ESource
)
{
    if (!extESource_.valid())
    {
        extESource_ =
            tmp<fvScalarMatrix>
            (
                new fvScalarMatrix(e(), dimEnergy/dimTime)
            );
    }
    extESource_.ref() += ESource;
}


void Foam::phaseCompressibleSystem::addUCoeff
(
    const volScalarField::Internal& UCoeff
)
{
    if (!dragSource_.valid())
    {
        dragSource_ = tmp<fvVectorMatrix>(new fvVectorMatrix(U_, dimForce));
    }
    dragSource_.ref() -= fvm::Sp(UCoeff, U_);
}


void Foam::phaseCompressibleSystem::addUSource(const volVectorField::Internal& USource)
{
    if (!dragSource_.valid())
    {
        dragSource_ = tmp<fvVectorMatrix>(new fvVectorMatrix(U_, dimForce));
    }
    dragSource_.ref() += USource;
}


const Foam::blast::momentumTransportModel&
Foam::phaseCompressibleSystem::turbulence() const
{
    return turbulence_();
}


Foam::blast::momentumTransportModel&
Foam::phaseCompressibleSystem::turbulence()
{
    return turbulence_();
}


Foam::scalar Foam::phaseCompressibleSystem::CoNum() const
{
    surfaceScalarField amaxSf
    (
        fvc::interpolate(speedOfSound())*mesh_.magSf()
    );
    // Remove wave speed from wedge boundaries
    forAll(amaxSf.boundaryField(), patchi)
    {
        if (isA<wedgeFvPatch>(mesh_.boundary()[patchi]))
        {
            amaxSf.boundaryFieldRef() = Zero;
        }
    }
    amaxSf += mag(fvc::flux(U_));

    scalarField sumAmaxSf
    (
        fvc::surfaceSum(amaxSf)().primitiveField()
    );

    return 0.5*gMax(sumAmaxSf/mesh_.V().field())*mesh_.time().deltaTValue();
}


bool Foam::phaseCompressibleSystem::writeData(Ostream& os) const
{
    return os.good();
}

bool Foam::phaseCompressibleSystem::read()
{
    return regIOobject::read();
}

// ************************************************************************* //
