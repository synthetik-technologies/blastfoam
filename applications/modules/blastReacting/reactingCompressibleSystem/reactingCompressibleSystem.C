/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2025
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

#include "reactingCompressibleSystem.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static member functions * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(reactingCompressibleSystem, 0);
    addToRunTimeSelectionTable
    (
        compressibleSystem,
        reactingCompressibleSystem,
        singlePhase
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactingCompressibleSystem::reactingCompressibleSystem
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    compressibleSystem(dict, mesh),
    thermo_(psiMulticomponentThermo::New(mesh)),
    rho_
    (
        IOobject
        (
            "rho",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rho", dimDensity, 1.0)
    ),
    p_(thermo_->p()),
    T_(thermo_->T()),
    e_(thermo_->he()),
    odeCombustion_(false)
{
    thermo_->validate("compressibleSystem", "e");
    rho_ = thermo_->rho();

    turbulence_ =
        compressible::momentumTransportModel::New
        (
            rho_,
            U_,
            rhoPhi_,
            thermo_()
        );

    mesh.schemes().setFluxRequired(U_.name());

    thermophysicalTransport_ =
        fluidMulticomponentThermophysicalTransportModel::New
        (
            turbulence_(),
            thermo_()
        );

    if (thermo_->Y().size() > 1)
    {
        reaction_.set
        (
            combustionModel::New
            (
                thermo_(),
                turbulence_()
            ).ptr()
        );
        odeCombustion_ = dict.lookupOrDefault("odeCombustion", false);
    }

    typeIOobject<IOdictionary> radIO
    (
        "radiationProperties",
        mesh.time().constant(),
        mesh
    );
    if (radIO.headerOk())
    {
        radiation_ = radiationModel::New(T_);
    }
    else
    {
        dictionary radDict;
        radDict.add("radiationModel", "none");
        radiation_ = radiationModel::New(radDict, T_);
    }

    fluxScheme_ = fluxScheme::NewSingle(phi_);
    encode();

    // Mark flux fields to be cached
    mesh.addTemporaryObject(reconstruction::ownName(rho_.name()));
    mesh.addTemporaryObject(reconstruction::neiName(rho_.name()));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reactingCompressibleSystem::~reactingCompressibleSystem()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::reactingCompressibleSystem::solve()
{
    if (reaction_.valid() && (odeCombustion_ || this->step() == 1))
    {
        reaction_->correct();
    }

    // Save current density for "old" value of specie masses
    tmp<volScalarField> trho0;
    if (reaction_.valid())
    {
        trho0 = volScalarField::New("rho0", rho_);
    }

    volScalarField deltaRho(fvc::div(rhoPhi_));
    volVectorField deltaRhoU(fvc::div(rhoUPhi_));
    volScalarField deltaRhoE(fvc::div(rhoEPhi_));
    if (reaction_.valid())
    {
        deltaRhoE -= reaction_->Qdot();
    }

    addSources(deltaRhoU, deltaRhoE);

    //- Store and blend old values
    this->storeAndBlendOld(rho_);
    this->storeAndBlendOld(rhoU_);
    this->storeAndBlendOld(rhoE_);

    //- Store changed in mass, momentum and energy
    this->storeAndBlendDelta(deltaRho);
    this->storeAndBlendDelta(deltaRhoU);
    this->storeAndBlendDelta(deltaRhoE);


    const dimensionedScalar& dT = rho_.time().deltaT();
    rho_ -= dT*deltaRho;
    rhoU_ -= dT*deltaRhoU;
    rhoE_ -= dT*deltaRhoE;

    if (reaction_.valid())
    {
        const volScalarField& rho0 = trho0();
        PtrList<volScalarField>& Ys = thermo_->Y();
        forAll(Ys, i)
        {
            if (thermo_->solveSpecie(i))
            {
                volScalarField deltaRhoY
                (
                    fvc::div(rhoPhi_, Ys[i], "div(" + rhoPhi_.name() + ",Yi)")
                );
                deltaRhoY.internalFieldRef() -= reaction_->R(i);

                volScalarField rhoYiOld(rho0*Ys[i]);
                this->storeAndBlendOld(rhoYiOld);
                this->storeAndBlendDelta(deltaRhoY);

                Ys[i] = (rhoYiOld - dT*deltaRhoY)/rho_;
                Ys[i].max(0.0);
                Ys[i].correctBoundaryConditions();
            }
        }
    }
}


void Foam::reactingCompressibleSystem::postImplicit()
{
    this->decode();

    // Solve mass
    if (needSolve(rho_.name()))
    {
        fvScalarMatrix rhoEqn
        (
            fvm::ddt(rho_) - rhoAdvection_()
        ==
            models().source(rho_)
        );

        constraints().constrain(rhoEqn);
        rhoEqn.solve();
        constraints().constrain(rho_);
    }

    if (reaction_.valid())
    {
        Info<< "Solving reactions" << endl;

        PtrList<volScalarField>& Y = thermo_->Y();
        forAll(Y, i)
        {
            if (thermo_->solveSpecie(i))
            {
                volScalarField& Yi = Y[i];
                fvScalarMatrix YiEqn
                (
                    fvm::ddt(rho_, Yi) - rhoYAdvection_[i]()
                  + thermophysicalTransport_->divj(Yi)
                 ==
                    // reaction_->R(Yi)
                    models().source(rho_, Yi)
                );

                constraints().constrain(YiEqn);
                YiEqn.solve("Yi");
                constraints().constrain(Yi);

                Yi.max(0.0);
            }
        }
        thermo_->normaliseY();
    }

    // Viscous terms if not using explicit viscosity
    compressibleSystem::postImplicit();

    // Update thermo
    thermo_->correct();

    p_.internalFieldRef() = rho_()/thermo_->psi()();
    constraints().constrain(p_);
    p_.correctBoundaryConditions();

    // Update density boundary conditions
    rho_.boundaryFieldRef() ==
        thermo_->psi().boundaryField()*p_.boundaryField();
}


void Foam::reactingCompressibleSystem::update()
{
    decode();
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


void Foam::reactingCompressibleSystem::decode()
{
    thermo_->normaliseY();

    thermo_->rho() = rho_;

    U_.internalFieldRef() = rhoU_()/rho_();
    U_.correctBoundaryConditions();

    K_ = 0.5*magSqr(U_);

    e_.internalFieldRef() = rhoE_()/rho_() - K_();
    e_.correctBoundaryConditions();

    thermo_->correct();
    p_.internalFieldRef() = rho_/thermo_->psi();
    p_.correctBoundaryConditions();
    rho_.boundaryFieldRef() ==
        thermo_->psi().boundaryField()*p_.boundaryField();

    rhoU_.boundaryFieldRef() = rho_.boundaryField()*U_.boundaryField();
    rhoE_.boundaryFieldRef() =
        rho_.boundaryField()*(e_.boundaryField() + K_.boundaryField());
}


void Foam::reactingCompressibleSystem::storeExplicit()
{
    compressibleSystem::storeExplicit();

    rhoAdvection_ = fvc::ddt(rho_);

    if (reaction_.valid())
    {
        PtrList<volScalarField>& Y = thermo_->Y();
        rhoYAdvection_.setSize(Y.size());
        forAll(Y, phasei)
        {
            rhoYAdvection_[phasei] = fvc::ddt(rho_, Y[phasei]);
        }
    }
}


void Foam::reactingCompressibleSystem::clear()
{
    compressibleSystem::clear();

    rhoAdvection_.clear();
    forAll(rhoYAdvection_, phasei)
    {
        rhoYAdvection_[phasei].clear();
    }
}


void Foam::reactingCompressibleSystem::encode()
{
    K_ = 0.5*magSqr(U_);

    rho_ = thermo_->rho();
    rhoU_ = rho_*U_;
    rhoE_ = rho_*(e_ + K_);
}


Foam::tmp<Foam::volScalarField>
Foam::reactingCompressibleSystem::speedOfSound() const
{
    return sqrt(thermo_->Cp()/(thermo_->Cv()*thermo_->psi()));
}

// ************************************************************************* //
