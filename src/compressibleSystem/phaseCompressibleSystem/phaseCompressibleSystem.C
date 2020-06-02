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
            blastCompressibleTurbulenceModel::New
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
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
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

void Foam::phaseCompressibleSystem::solve
(
    const label stepi,
    const scalarList& ai,
    const scalarList& bi
)
{
    const fvMesh& mesh(rho_.mesh());
    volVectorField rhoUOld(rhoU_);
    volScalarField rhoEOld(rhoE_);
    if (mesh.moving() && stepi == 1)
    {
        volScalarField::Internal v0Byv(mesh.Vsc0()/mesh.Vsc());
        rhoUOld.ref() *= v0Byv;
        rhoEOld.ref() *= v0Byv;
    }

    if (oldIs_[stepi - 1] != -1)
    {
        rhoUOld_.set
        (
            oldIs_[stepi - 1],
            new volVectorField(rhoUOld)
        );
        rhoEOld_.set
        (
            oldIs_[stepi - 1],
            new volScalarField(rhoEOld)
        );


    }
    rhoUOld *= ai[stepi - 1];
    rhoEOld *= ai[stepi - 1];

    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = oldIs_[i];
        if (fi != -1 && ai[fi] != 0)
        {
            rhoUOld += ai[fi]*rhoUOld_[fi];
            rhoEOld += ai[fi]*rhoEOld_[fi];
        }
    }

    volVectorField deltaRhoU(fvc::div(rhoUPhi_) - g_*rho_);
    volScalarField deltaRhoE
    (
        fvc::div(rhoEPhi_)
      - ESource()
      - (rhoU_ & g_)
    );
    if (extESource_.valid())
    {
        deltaRhoE.ref() += extESource_();
    }

    if (deltaIs_[stepi - 1] != -1)
    {
        deltaRhoU_.set
        (
            deltaIs_[stepi - 1],
            new volVectorField(deltaRhoU)
        );
        deltaRhoE_.set
        (
            deltaIs_[stepi - 1],
            new volScalarField(deltaRhoE)
        );
    }

    scalar f(bi[stepi - 1]);
    deltaRhoU *= bi[stepi - 1];
    deltaRhoE *= bi[stepi - 1];

    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = deltaIs_[i];
        if (fi != -1 && bi[fi] != 0)
        {
            f += bi[fi];
            deltaRhoU += bi[fi]*deltaRhoU_[fi];
            deltaRhoE += bi[fi]*deltaRhoE_[fi];
        }
    }

    dimensionedScalar dT = rho_.time().deltaT();
    vector solutionDs((vector(rho_.mesh().solutionD()) + vector::one)/2.0);
    rhoU_ = cmptMultiply(rhoUOld - dT*deltaRhoU, solutionDs);
    rhoE_ = rhoEOld - dT*deltaRhoE;
    if (radiation_->type() != "none")
    {
        calcAlphaAndRho();
        e() = rhoE_/rho_ - 0.5*magSqr(U_);
        e().correctBoundaryConditions();
        rhoE_ =
            radiation_->calcRhoE
            (
                f*dT,
                rhoE_,
                rho_,
                e(),
                this->thermo().Cv()
            );
    }

    if (stepi == oldIs_.size())
    {
        radiation_->correct();
    }

    if
    (
        stepi == oldIs_.size()
     && (
            turbulence_.valid()
         || dragSource_.valid()
        )
    )
    {
        calcAlphaAndRho();
        U_ = rhoU_/rho_;
        U_.correctBoundaryConditions();

        fvVectorMatrix UEqn
        (
            fvm::ddt(rho_, U_) - fvc::ddt(rho_, U_)
        );
        if (turbulence_.valid())
        {
            volScalarField muEff("muEff", turbulence_->muEff());
            volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U_))));
            UEqn -=
                fvm::laplacian(muEff, U_)
              + fvc::div(tauMC);
        }
        if (dragSource_.valid())
        {
            UEqn += dragSource_();
        }
        UEqn.solve();

        rhoU_ = rho_*U_;

        if (turbulence_.valid())
        {
            e() = rhoE_/rho_ - 0.5*magSqr(U_);
            e().correctBoundaryConditions();

            Foam::solve
            (
                fvm::ddt(rho_, e()) - fvc::ddt(rho_, e())
            - fvm::laplacian(turbulence_->alphaEff(), e())
            );
            rhoE_ = rho_*(e() + 0.5*magSqr(U_)); // Includes change to total energy from viscous term in momentum equation

            turbulence_->correct();
        }
    }
}


void Foam::phaseCompressibleSystem::setODEFields
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
    rhoUOld_.resize(nOld_);
    rhoEOld_.resize(nOld_);

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
    deltaRhoU_.resize(nDelta_);
    deltaRhoE_.resize(nDelta_);
}

void Foam::phaseCompressibleSystem::clearODEFields()
{
    fluxScheme_->clear();
    rhoUOld_.clear();
    rhoEOld_.clear();
    rhoUOld_.resize(nOld_);
    rhoEOld_.resize(nOld_);

    deltaRhoU_.clear();
    deltaRhoE_.clear();
    deltaRhoU_.resize(nDelta_);
    deltaRhoE_.resize(nDelta_);

    extESource_.clear();
    dragSource_.clear();
}


void Foam::phaseCompressibleSystem::addESource(const volScalarField::Internal& extEsrc)
{
    if (!extESource_.valid())
    {
        extESource_ =
            tmp<volScalarField::Internal>(new volScalarField::Internal("extEsrc", extEsrc));
    }
    else
    {
        extESource_.ref() += extEsrc;
    }
}


void Foam::phaseCompressibleSystem::addUCoeff(const volScalarField::Internal& UCoeff)
{
    if (!dragSource_.valid())
    {
        dragSource_ = tmp<fvVectorMatrix>(new fvVectorMatrix(U_, dimForce));
    }
    dragSource_.ref() += fvm::Sp(UCoeff, U_);
}


void Foam::phaseCompressibleSystem::addUSource(const volVectorField::Internal& USource)
{
    if (!dragSource_.valid())
    {
        dragSource_ = tmp<fvVectorMatrix>(new fvVectorMatrix(U_, dimForce));
    }
    dragSource_.ref() += USource;
}


const Foam::blastCompressibleTurbulenceModel&
Foam::phaseCompressibleSystem::turbulence() const
{
    return turbulence_();
}


Foam::blastCompressibleTurbulenceModel&
Foam::phaseCompressibleSystem::turbulence()
{
    return turbulence_();
}


bool Foam::phaseCompressibleSystem::writeData(Ostream& os) const
{
    return os.good();
}

// ************************************************************************* //
