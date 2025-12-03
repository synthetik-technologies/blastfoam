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

#include "compressibleSystem.H"
#include "uniformDimensionedFields.H"
#include "fvm.H"
#include "fvcSmooth.H"
#include "MULES.H"
#include "fvcMeshPhi.H"
#include "wedgeFvPatch.H"
#include "emptyFvPatch.H"
#include "fluidThermoThermophysicalTransportModel.H"
#include "fluidMulticomponentThermophysicalTransportModel.H"
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
        {
            dictionary turbulenceDict;
            {
                typeIOobject<IOdictionary> momentumTransport
                (
                    momentumTransportModel::typeName,
                    this->mesh().time().constant(),
                    this->mesh(),
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE,
                    false
                );

                if (momentumTransport.headerOk())
                {
                    turbulenceDict = IOdictionary(momentumTransport);
                }
                else
                {
                    typeIOobject<IOdictionary> turbulenceProperties
                    (
                        "turbulenceProperties",
                        this->mesh().time().constant(),
                        this->mesh(),
                        IOobject::MUST_READ_IF_MODIFIED,
                        IOobject::NO_WRITE,
                        false
                    );

                    if (turbulenceProperties.headerOk())
                    {
                        turbulenceDict = IOdictionary(turbulenceProperties);
                    }
                    else
                    {
                        turbulenceDict = IOdictionary(momentumTransport);
                    }
                }
            }
            const word modelType
            (
                turbulenceDict.lookup("simulationType")
            );

            Info<< "Selecting turbulence model type " << modelType << endl;

            compressibleMomentumTransportModel::
            dictionaryConstructorTable::iterator cstrIter =
                compressibleMomentumTransportModel::
                dictionaryConstructorTablePtr_->find(modelType);

            if
            (
                cstrIter
             == compressibleMomentumTransportModel::
                dictionaryConstructorTablePtr_->end()
            )
            {
                FatalErrorInFunction
                    << "Unknown "
                    << compressibleMomentumTransportModel::typeName
                    << " type "
                    << modelType << nl << nl
                    << "Valid "
                    << compressibleMomentumTransportModel::typeName
                    << " types:" << endl
                    << compressibleMomentumTransportModel::
                        dictionaryConstructorTablePtr_->sortedToc()
                    << exit(FatalError);
            }

            turbulence_ = autoPtr<compressibleMomentumTransportModel>
            (
                cstrIter()
                (
                    geometricOneField(),
                    rhoEff(),
                    U(),
                    rhoPhi(),
                    phi(),
                    thermo()
                )
            );
        }
        // turbulence_ =
        //     momentumTransportModel::New<compressibleMomentumTransportModel>
        //     (
        //         geometricOneField(),
        //         rhoEff(),
        //         U(),
        //         rhoPhi(),
        //         phi(),
        //         this->thermo()
        //     );
        turbulence_->validate();

        mesh().schemes().setFluxRequired(U_.name());

        if (isA<multicomponentThermo>(this->thermo()))
        {
            thermophysicalTransport_ =
            (
                fluidMulticomponentThermophysicalTransportModel::New
                (
                    turbulence_,
                    dynamicCast<const fluidMulticomponentThermo>
                    (
                        this->thermo()
                    )
                ).ptr()
            );
        }
        else
        {
            thermophysicalTransport_ =
            (
                fluidThermoThermophysicalTransportModel::New
                (
                    turbulence_,
                    this->thermo()
                ).ptr()
            );
        }
    }
    else
    {
        explicitViscosity_ = false;
    }
}


void Foam::compressibleSystem::addSources
(
    volVectorField::Internal& rhoUSource,
    volScalarField::Internal& rhoESource
) const
{

    if (mag(g_).value() > small)
    {
        rhoUSource -= g_*rhoEff()();
        rhoESource -= g_ & rhoU_();
    }

    if (explicitViscosity_ && turbulence_.valid())
    {
        tmp<volSymmTensorField> tdevTau(turbulence_->devTau());
        // tdevTau.ref() += (2.0/3.0)*rhoEff()*turbulence_->k()*symmTensor::I;

        rhoUSource += fvc::div(tdevTau());
        rhoESource +=
            fvc::div
            (
                fvc::dotInterpolate(mesh().Sf(), tdevTau)
              & flux().Uf()
            )
          + fvc::div(thermophysicalTransport_->q()*mesh().magSf());
    }
}



void Foam::compressibleSystem::updateCorDeltaT()
{
    if (fv::localEulerDdt::enabled(mesh()))
    {
        if (!corDeltaTPtr_.valid())
        {
            corDeltaTPtr_.set
            (
                new volScalarField
                (
                    IOobject
                    (
                        "corDeltaT",
                        mesh().time().name(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh(),
                    1.0
                )
            );
            localRDeltaTPtr_.set
            (
                new volScalarField
                (
                    IOobject
                    (
                        fv::localEulerDdt::rDeltaTName,
                        mesh().time().name(),
                        mesh()
                    ),
                    mesh(),
                    1.0/mesh().time().deltaT()
                )
            );
        }

        const dimensionedScalar& deltaT = mesh().time().deltaT();
        const surfaceScalarField& magSf = mesh().magSf();
        surfaceScalarField amaxSf(fvc::interpolate(speedOfSound())*magSf);

        // Remove wave speed from wedge boundaries
        surfaceScalarField::Boundary& bamaxSf = amaxSf.boundaryFieldRef();
        forAll(bamaxSf, patchi)
        {
            if (isA<wedgeFvPatch>(mesh().boundary()[patchi]))
            {
                bamaxSf[patchi] = Zero;
            }
        }
        amaxSf += mag(this->phi());

        const dictionary& pimpleDict =
                mesh().solution().subOrEmptyDict("PIMPLE");
        const scalar maxCo =
            pimpleDict.lookupOrDefault("maxCo", fvTimeInt_->maxCo()/2.0);

        // Calulate rDeltaT (local)
        volScalarField& rDeltaT = localRDeltaTPtr_();
        rDeltaT.internalFieldRef() =
            fvc::surfaceSum(amaxSf)()()/((2*maxCo)*mesh().V());

        if (explicitViscosity_)
        {
            surfaceScalarField deltaCoeffSqr(magSf*mesh().deltaCoeffs());

            // Remove wave speed from wedge boundaries
            surfaceScalarField::Boundary& bdeltaCoeffSqr =
                deltaCoeffSqr.boundaryFieldRef();
            forAll(bdeltaCoeffSqr, patchi)
            {
                if (isA<wedgeFvPatch>(mesh().boundary()[patchi]))
                {
                    bdeltaCoeffSqr[patchi] = Zero;
                }
            }

            rDeltaT.internalFieldRef() =
            (
                max
                (
                    rDeltaT.internalField(),
                    fvc::surfaceSum
                    (
                        deltaCoeffSqr
                       *fvc::interpolate(turbulence_->nuEff())
                    )()()
                   /(mesh().V())
                )
            );
        }

        scalar minRDeltaT(gMin(rDeltaT.primitiveField()));
        if (pimpleDict.found("maxDeltaT"))
        {
            const scalar clipRDeltaT =
                1.0/pimpleDict.lookup<scalar>("maxDeltaT");
            rDeltaT.max(clipRDeltaT);
            minRDeltaT = max(minRDeltaT, clipRDeltaT);
        }
        scalar maxRDeltaT(gMax(rDeltaT.primitiveField()));
        if (pimpleDict.found("minDeltaT"))
        {
            const scalar clipRDeltaT =
                1.0/pimpleDict.lookup<scalar>("minDeltaT");
            rDeltaT.min(clipRDeltaT);
            maxRDeltaT = min(maxRDeltaT, clipRDeltaT);
        }
        rDeltaT.correctBoundaryConditions();

        Info<< "Flow time scale min/max = "
            << 1.0/maxRDeltaT << ", "
            << 1.0/minRDeltaT << endl;


        const scalar rDeltaTSmoothingCoeff =
            pimpleDict.lookupOrDefault("rDeltaTSmoothingCoeff", 0.02);
        if (rDeltaTSmoothingCoeff > 0)
        {
            fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);

            Info<< "Smoothed flow time scale min/max = "
                << 1.0/gMax(rDeltaT.primitiveField()) << ", "
                << 1.0/gMin(rDeltaT.primitiveField()) << endl;
        }

        volScalarField& corDeltaT = corDeltaTPtr_();
        corDeltaT = rDeltaT*deltaT;
    }
}


void Foam::compressibleSystem::addUSource(fvVectorMatrix& UEqn) const
{
    if (dragSource_.valid())
    {
        UEqn -= dragSource_;
    }
}

void Foam::compressibleSystem::addESource(fvScalarMatrix& EEqn) const
{
    if (dragSource_.valid())
    {
        EEqn -= (dragSource_ & U_) & U_;
    }
    if (extESource_.valid())
    {
        EEqn -= extESource_;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressibleSystem::compressibleSystem
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    timeIntegrationSystem("compressibleSystem", mesh),
    U_
    (
        IOobject
        (
            "U",
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    K_
    (
        IOobject
        (
            "K",
            mesh.time().name(),
            mesh
        ),
        0.5*magSqr(U_)
    ),
    rhoU_
    (
        IOobject
        (
            "rhoU",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimDensity*dimVelocity, Zero)
    ),
    rhoE_
    (
        IOobject
        (
            "rhoE",
            mesh.time().name(),
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
            mesh.time().name(),
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
            mesh.time().name(),
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
            mesh.time().name(),
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
            mesh.time().name(),
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
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        ),
        dimensionedVector(dimAcceleration, Zero)
    ),
    solutionDs_((vector(mesh.solutionD()) + vector::one)/2.0),
    explicitViscosity_(dict.lookupOrDefault("explicitViscosity", false))
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
    K_ = 0.5*magSqr(U_);
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
    if (step() == 0) updateCorDeltaT();
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


void Foam::compressibleSystem::postExplicit()
{}

void Foam::compressibleSystem::postImplicit()
{
    if (turbulence_.valid())
    {
        turbulence_->predict();
    }
    if (thermophysicalTransport_.valid())
    {
        thermophysicalTransport_->predict();
    }

    tmp<surfaceVectorField> devTau;
    if (needSolve_U())
    {
        tmp<fvVectorMatrix> divDevTau;
        if (!explicitViscosity_ && turbulence_.valid())
        {
            divDevTau =
                turbulence_->divDevTau(U_);
                // + fvc::grad((2.0/3.0)*rhoEff()*turbulence_->k());
        }

        // Solve momentum
        fvVectorMatrix UEqn
        (
            fvm::ddt(rhoEff(), U_) - rhoUAdvection_()
        ==
            models().source(rhoEff(), U_)
        );

        if (divDevTau.valid())
        {
            UEqn += divDevTau();
        }
        addUSource(UEqn);

        UEqn.relax();

        constraints().constrain(UEqn);
        UEqn.solve();
        constraints().constrain(U_);

        if (divDevTau.valid())
        {
            devTau = divDevTau().flux();
        }

        // Update kinetic energy and momentum
        K_ = 0.5*magSqr(U_);
        rhoU_ = rhoEff()*U_;
    }

    // Solve thermal energy diffusion
    if (needSolve_E())
    {
        volScalarField& he = thermo().he();
        fvScalarMatrix EEqn
        (
            fvm::ddt(rhoEff(), he)
          - rhoEAdvection_()        // Explicit advection contribtion
          + fvc::ddt(rhoEff(), K_)  // Change in kinetic energy
         ==
            models().source(rhoEff(), he)
        );

        if (devTau.valid())
        {
            EEqn +=
                fvc::div(devTau & flux().Uf())
              + thermophysicalTransport_->divq(he);
        }

        EEqn.relax();

        constraints().constrain(EEqn);
        EEqn.solve();
        constraints().constrain(he);

        // Update total energy
        rhoE_ = rhoEff()*(he + K_);
    }

    if (turbulence_.valid())
    {
        turbulence_->correct();
    }

    if (thermophysicalTransport_.valid())
    {
        thermophysicalTransport_->correct();
    }
}


Foam::volScalarField& Foam::compressibleSystem::rhoEff()
{
    return rho();
}


const Foam::volScalarField& Foam::compressibleSystem::rhoEff() const
{
    return rho();
}


void Foam::compressibleSystem::storeExplicit()
{
    if (needSolve_U())
    {
        rhoUAdvection_ = fvc::ddt(rhoU_);
    }

    if (needSolve_E())
    {
        rhoEAdvection_ = fvc::ddt(rhoE_);
    }
}


void Foam::compressibleSystem::clear()
{
    fluxScheme_->clear();
    rhoEAdvection_.clear();
    rhoUAdvection_.clear();
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
    // if (mesh().topoChanged())
    // {
    //     const_cast<compressibleSystem&>(*this).decode();
    // }
    const surfaceScalarField& magSf = mesh().magSf();
    surfaceScalarField amaxSf(fvc::interpolate(speedOfSound())*magSf);

    // Remove wave speed from wedge boundaries
    surfaceScalarField::Boundary& bamaxSf = amaxSf.boundaryFieldRef();
    forAll(bamaxSf, patchi)
    {
        if (isA<wedgeFvPatch>(mesh().boundary()[patchi]))
        {
            bamaxSf[patchi] = Zero;
        }
    }
    amaxSf += mag(fvc::flux(this->U()));

    scalarField sumAmaxSf
    (
        fvc::surfaceSum(amaxSf)().primitiveField()
    );

    tmp<volScalarField::Internal> tV(mesh().Vsc());
    const scalarField& V = tV();

    scalarField cof(0.5*(sumAmaxSf/V)*mesh().time().deltaTValue());
    scalar CoNum = 0.5*gMax(sumAmaxSf/V)*mesh().time().deltaTValue();

    scalar meanCoNum =
        0.5
       *(gSum(sumAmaxSf)/gSum(V))
       *mesh().time().deltaTValue();

    Info<< "Courant Number ";
    if (mesh().name() != polyMesh::defaultRegion)
    {
        Info<< "for region " << mesh().name() << " ";
    }
    Info<< "Mean = " << meanCoNum << ", Max = "<< CoNum << endl;

    return CoNum;
}


Foam::scalar Foam::compressibleSystem::DiNum() const
{
    // Check diffusion number
    scalar DiNum = 0.0;
    if (explicitViscosity_)
    {
        const surfaceScalarField& magSf = mesh().magSf();
        surfaceScalarField deltaCoeffSqr(magSf*mesh().deltaCoeffs());

        // Remove wave speed from wedge boundaries
        surfaceScalarField::Boundary& bdeltaCoeffSqr =
            deltaCoeffSqr.boundaryFieldRef();
        forAll(bdeltaCoeffSqr, patchi)
        {
            if (isA<wedgeFvPatch>(mesh().boundary()[patchi]))
            {
                bdeltaCoeffSqr[patchi] = Zero;
            }
        }

        const volScalarField::Internal DiNumvf
        (
            fvc::surfaceSum
            (
                deltaCoeffSqr
               *fvc::interpolate(turbulence_->nuEff())
            )()()
           /(mesh().V())
           *mesh().time().deltaT()
        );
        const scalar meanDiNum = gAverage(DiNumvf);
        const scalar maxDiNum = gMax(DiNumvf);

        Info<< "Diffusion Number mean: " << meanDiNum
            << " max: " << maxDiNum << endl;

        DiNum = maxDiNum;
    }
    return DiNum;
}

// ************************************************************************* //
