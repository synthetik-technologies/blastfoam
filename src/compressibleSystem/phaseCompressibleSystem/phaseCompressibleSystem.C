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
#include "blastCompressibleTurbulenceModel.H"
#include "uniformDimensionedFields.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseCompressibleSystem, 0);
    defineRunTimeSelectionTable(phaseCompressibleSystem, dictionary);
}


void Foam::phaseCompressibleSystem::setModels(const dictionary& dict)
{
    if (Foam::max(this->thermo().mu()).value() > 0)
    {
        turbulence_ =
        (
            blast::turbulenceModel::New
            (
                rho_,
                U_,
                rhoPhi_,
                this->thermo()
            )
        );
        turbulence_->validate();
    }

    bool usesRadProperties = false;
    {
        IOdictionary radIODict
        (
            IOobject
            (
                "radiationProperties",
                rho_.time().constant(),
                rho_.mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        );
        usesRadProperties = radIODict.found("radiationModel");
    }


    if (dict.found("radiationModel"))
    {
        radiation_.set(radiationModel::New(dict, this->T()).ptr());
    }
    else if(usesRadProperties)
    {
        radiation_.set(radiationModel::New(this->T()).ptr());
    }
    else
    {
        dictionary radDict;
        radDict.add("radiationModel", "none");
        radiation_.set(radiationModel::New(radDict, this->T()).ptr());
    }

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseCompressibleSystem::phaseCompressibleSystem
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
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
        dimensionedScalar("rho", dimDensity, 0.0),
        wordList(mesh.boundaryMesh().size(), "zeroGradient")
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
        wordList(p_.boundaryField().types().size(), "zeroGradient")
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
        dimensionedScalar("0", dimDensity*sqr(dimVelocity), 0.0),
        wordList(p_.boundaryField().types().size(), "zeroGradient")
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
    TLow_("TLow", dimTemperature, 0.0)
{
    TLow_.readIfPresent(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseCompressibleSystem::~phaseCompressibleSystem()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::phaseCompressibleSystem::solve()
{
    if (step() == 1)
    {
        rhoOldTmp_ = tmp<volScalarField>(new volScalarField(rho_));
    }

    volVectorField rhoUOld(rhoU_);
    volScalarField rhoEOld(rhoE_);

    //- Store old values
    this->storeAndBlendOld(rhoUOld, rhoUOld_);
    this->storeAndBlendOld(rhoEOld, rhoEOld_);

    //- Calculate deltas for momentum and energy
    volVectorField deltaRhoU
    (
        fvc::div(rhoUPhi_) - g_*rho_
    );

    volScalarField deltaRhoE
    (
        fvc::div(rhoEPhi_)
      - ESource()
      - (rhoU_ & g_)
    );

    //- Store changed in momentum and energy
    this->storeAndBlendDelta(deltaRhoU, deltaRhoU_);
    this->storeAndBlendDelta(deltaRhoE, deltaRhoE_);

    //- Solve for momentum and energy
    dimensionedScalar dT = rho_.time().deltaT();
    vector solutionDs((vector(rho_.mesh().solutionD()) + vector::one)/2.0);
    rhoU_ = cmptMultiply(rhoUOld - dT*deltaRhoU, solutionDs);
    rhoE_ = rhoEOld - dT*deltaRhoE;
    Info<< "energy: "
        << sum(rhoE_()*rho_.mesh().V()).value()<<endl;
}


void Foam::phaseCompressibleSystem::postUpdate()
{
    this->decode();

    if (radiation_->type() != "none")
    {
        radiation_->correct();

        calcAlphaAndRho();
        U_ = rhoU_/rho_;
        U_.correctBoundaryConditions();

        e() = rhoE_/rho_ - 0.5*magSqr(U_);
        e().correctBoundaryConditions();
        rhoE_ =
            radiation_->calcRhoE
            (
                rho_.mesh().time().deltaT(),
                rhoE_,
                rho_,
                e(),
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
        rho_.oldTime() = rhoOldTmp_();

        fvVectorMatrix UEqn
        (
            fvm::ddt(rho_, U_) - fvc::ddt(rho_, U_)
        );
        fvScalarMatrix eEqn
        (
            fvm::ddt(rho_, e()) - fvc::ddt(rho_, e())
        );

        if (dragSource_.valid())
        {
            UEqn -= dragSource_();
        }
        if (extESource_.valid())
        {
            eEqn -= extESource_();
        }


        if (turbulence_.valid())
        {
            UEqn += turbulence_->divDevRhoReff(U_);
            eEqn -= fvm::laplacian(turbulence_->alphaEff(), e());
        }

        UEqn.solve();
        eEqn.solve();

        rhoU_ = rho_*U_;
        rhoE_ = rho_*(e() + 0.5*magSqr(U_)); // Includes change to total energy from viscous term in momentum equation
    }

    if (turbulence_.valid())
    {
        turbulence_->correct();
    }

    this->thermo().postUpdate();
}


void Foam::phaseCompressibleSystem::clearODEFields()
{
    fluxScheme_->clear();

    this->clearOld(rhoUOld_);
    this->clearOld(rhoEOld_);

    this->clearDelta(deltaRhoU_);
    this->clearDelta(deltaRhoE_);

    extESource_.clear();
    dragSource_.clear();

    if (rhoOldTmp_.valid())
    {
        rhoOldTmp_.clear();
    }
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


const Foam::blast::turbulenceModel&
Foam::phaseCompressibleSystem::turbulence() const
{
    return turbulence_();
}


Foam::blast::turbulenceModel&
Foam::phaseCompressibleSystem::turbulence()
{
    return turbulence_();
}


bool Foam::phaseCompressibleSystem::writeData(Ostream& os) const
{
    return os.good();
}

// ************************************************************************* //
