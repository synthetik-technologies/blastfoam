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
#include "MULES.H"
#include "fvcMeshPhi.H"
#include "wedgeFvPatch.H"
#include "blastRadiationModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(compressibleSystem, 0);
    defineRunTimeSelectionTable(compressibleSystem, singlePhase);
    defineRunTimeSelectionTable(compressibleSystem, twoPhase);
    defineRunTimeSelectionTable(compressibleSystem, multiphase);
    defineRunTimeSelectionTable(compressibleSystem, coupled);
}


// * * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * //

void Foam::compressibleSystem::setModels()
{
    if (Foam::max(this->thermo().mu()).value() > small)
    {
        turbulence_ =
        (
            compressible::momentumTransportModel::New
            (
                rhoEff(),
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
}


Foam::tmp<Foam::volVectorField> Foam::compressibleSystem::rhoUSource() const
{
    return g_*rhoEff();
}


Foam::tmp<Foam::volScalarField> Foam::compressibleSystem::rhoESource() const
{
    return rhoU_ & g_;
}

void Foam::compressibleSystem::limitAlphaRhoPhis
(
    const UPtrList<volScalarField>& alphas,
    const UPtrList<volScalarField>& rhos,
    UPtrList<surfaceScalarField>& alphaRhoPhis,
    const surfaceScalarField& phi,
    const surfaceScalarField& rhoPhi
)
{
    PtrList<surfaceScalarField> alphaRhoPhiUDs(alphas.size());
    // forAll(alphaRhoPhiUDs, phasei)
    // {
    //     alphaRhoPhiUDs.set
    //     (
    //         phasei,
    //         upwind<scalar>(mesh(), phi).flux(alphas[phasei]*rhos[phasei])
    //     );
    //     alphaRhoPhis[phasei] -= alphaRhoPhiUDs[phasei];
    // }

    {
        UPtrList<scalarField> alphaRhoPhisInternal(alphas.size());
        forAll(alphaRhoPhisInternal, phasei)
        {
            alphaRhoPhisInternal.set(phasei, &alphaRhoPhis[phasei]);
        }
        MULES::limitSum(alphaRhoPhisInternal);
    }

    const surfaceScalarField::Boundary& phibf = phi_.boundaryField();
    forAll(phibf, patchi)
    {
        if (phibf[patchi].coupled())
        {
            UPtrList<scalarField> alphaRhoPhisPatch(alphas.size());
            forAll(alphaRhoPhisPatch, phasei)
            {
                alphaRhoPhisPatch.set
                (
                    phasei,
                    &alphaRhoPhis[phasei].boundaryFieldRef()[patchi]
                );
            }
            MULES::limitSum(alphaRhoPhisPatch);
        }
    }

    // forAll(alphaRhoPhis, phasei)
    // {
    //     alphaRhoPhis[phasei] += alphaRhoPhiUDs[0];
    // }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressibleSystem::compressibleSystem
(
    const fvMesh& mesh
)
:
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
            IOobject::READ_IF_PRESENT,
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
            IOobject::READ_IF_PRESENT,
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
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
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
    g_
    (
        IOobject
        (
            "g",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        ),
        dimensionedVector(dimAcceleration, Zero)
    ),
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::compressibleSystem::~compressibleSystem()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::compressibleSystem::encode()
{
    rhoU_ = rhoEff()*U_;
    rhoE_ = rhoEff()*(he() + 0.5*magSqr(U_));
}


void Foam::compressibleSystem::update()
{
    decode();
    fluxScheme_->update
    (
        rhoEff(),
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
        fvc::div(rhoUPhi_) - g_*rhoEff()
    );
    this->fvTimeInt_->addDeltaSource(rhoU_.name(), deltaRhoU);

    volScalarField deltaRhoE
    (
        "deltaRhoE",
        fvc::div(rhoEPhi_)
      - (rhoU_ & g_)
    );
    this->fvTimeInt_->addDeltaSource(rhoE_.name(), deltaRhoE);

    //- Store old values
    this->storeAndBlendOld(rhoU_);
    this->storeAndBlendOld(rhoE_);

    //- Store changed in momentum and energy
    this->storeAndBlendDelta(deltaRhoU);
    this->storeAndBlendDelta(deltaRhoE);

    //- Solve for momentum and energy
    dimensionedScalar dT = rhoEff().time().deltaT();
    rhoU_ -= cmptMultiply(dT*deltaRhoU, solutionDs_);
    rhoE_ -= dT*deltaRhoE;
}


void Foam::compressibleSystem::postUpdate()
{
    // Solve momentum
    if (needSolve(U_.name()) || turbulence_.valid())
    {
        fvVectorMatrix UEqn
        (
            fvm::ddt(rhoEff(), U_) - fvc::ddt(rhoU_)
        ==
            models().source(rhoEff(), U_)
        );
        if (turbulence_.valid())
        {
            UEqn += turbulence_->divDevTau(U_);
            rhoE_ +=
                rhoEff().mesh().time().deltaT()
                *fvc::div
                (
                    fvc::dotInterpolate(rhoEff().mesh().Sf(), turbulence_->devTau())
                  & fluxScheme_->Uf()
                );
        }
        constraints().constrain(UEqn);
        UEqn.solve();
        constraints().constrain(U_);

        rhoU_ = rhoEff()*U_;

        //- Update internal energy
        he() = rhoE_/rhoEff() - 0.5*magSqr(U_);
    }

    // Solve thermal energy diffusion
    if (needSolve(he().name()) || turbulence_.valid())
    {
        fvScalarMatrix eEqn
        (
            fvm::ddt(rhoEff(), he()) - fvc::ddt(rhoEff().prevIter(), he())
        ==
            models().source(rhoEff(), he())
        );
        if (turbulence_.valid())
        {
            eEqn += thermophysicalTransport_->divq(he());
        }
        constraints().constrain(eEqn);
        eEqn.solve();
        constraints().constrain(he());

        rhoE_ = rhoEff()*(he() + 0.5*magSqr(U_));
    }

    if (turbulence_.valid())
    {
        turbulence_->correct();
    }
    this->thermo().correct();
    constraints().constrain(thermo().p());
    thermo().p().correctBoundaryConditions();
}


Foam::volScalarField& Foam::compressibleSystem::rhoEff()
{
    return rho();
}


const Foam::volScalarField& Foam::compressibleSystem::rhoEff() const
{
    return rho();
}


void Foam::compressibleSystem::clear()
{
    fluxScheme_->clear();
}


void Foam::compressibleSystem::addRhoCoeff
(
    const volScalarField::Internal& coeff
)
{
    if (!rhoSource_.valid())
    {
        rhoSource_ =
            tmp<fvScalarMatrix>
            (
                new fvScalarMatrix(this->rhoEff(), dimMass/dimTime)
            );
    }
    rhoSource_.ref() -= fvm::Sp(coeff, this->rhoEff());
}


void Foam::compressibleSystem::addRhoSource
(
    const volScalarField::Internal& src
)
{
    if (!rhoSource_.valid())
    {
        rhoSource_ =
            tmp<fvScalarMatrix>
            (
                new fvScalarMatrix(this->rhoEff(), dimMass/dimTime)
            );
    }
    rhoSource_.ref() += src;
}


void Foam::compressibleSystem::addUCoeff
(
    const volScalarField::Internal& coeff
)
{
    if (!dragSource_.valid())
    {
        dragSource_ = tmp<fvVectorMatrix>(new fvVectorMatrix(U_, dimForce));
    }
    dragSource_.ref() -= fvm::Sp(coeff, U_);
}


void Foam::compressibleSystem::addUSource
(
    const volVectorField::Internal& src
)
{
    if (!dragSource_.valid())
    {
        dragSource_ = tmp<fvVectorMatrix>(new fvVectorMatrix(U_, dimForce));
    }
    dragSource_.ref() += src;
}


void Foam::compressibleSystem::addECoeff
(
    const volScalarField::Internal& coeff
)
{
    if (!extESource_.valid())
    {
        extESource_ =
            tmp<fvScalarMatrix>
            (
                new fvScalarMatrix(this->he(), dimEnergy/dimTime)
            );
    }
    extESource_.ref() -= fvm::Sp(coeff, this->he());
}


void Foam::compressibleSystem::addESource
(
    const volScalarField::Internal& src
)
{
    if (!extESource_.valid())
    {
        extESource_ =
            tmp<fvScalarMatrix>
            (
                new fvScalarMatrix(this->he(), dimEnergy/dimTime)
            );
    }
    extESource_.ref() += src;
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
            amaxSf.boundaryFieldRef()[patchi] = Zero;
        }
    }
    amaxSf += mag(fvc::relative(fvc::flux(U()), U()));

    scalarField sumAmaxSf
    (
        fvc::surfaceSum(amaxSf)().primitiveField()
    );

    tmp<volScalarField::Internal> tV(mesh().Vsc());
    const volScalarField::Internal& V = tV();

    scalarField cof(0.5*(sumAmaxSf/V.field())*mesh().time().deltaTValue());
    scalar CoNum =
        0.5*gMax(sumAmaxSf/V.field())*mesh().time().deltaTValue();

    scalar meanCoNum =
        0.5*(gSum(sumAmaxSf)/gSum(V.field()))*mesh().time().deltaTValue();

    Info<< "Courant Number ";
    if (mesh().name() != polyMesh::defaultRegion)
    {
        Info<< "for region " << mesh().name() << " ";
    }
    Info<< "Mean = " << meanCoNum << ", Max = "<< CoNum << endl;
    return CoNum;
}

// ************************************************************************* //
