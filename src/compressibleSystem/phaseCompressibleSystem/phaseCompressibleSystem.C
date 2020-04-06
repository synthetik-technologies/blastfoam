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
#include "fiveEqnCompressibleTurbulenceModel.H"
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
    if (Foam::max(mu()).value() > 0)
    {
        turbulence_ =
        (
            fiveEqnCompressibleTurbulenceModel::New
            (
                rho_,
                U_,
                rhoPhi_,
                *this
            )
        );
        turbulence_->validate();
    }

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
    dictionary radDict;
    if (dict.found("radiationModel"))
    {
        radDict = dict;
    }
    else if(radIODict.found("radiationModel"))
    {
        radDict = radIODict;
    }
    else
    {
        radDict.add("radiationModel", "none");
    }
    radiation_.set(radiationModel::New(radDict, this->T()).ptr());
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
    fluxScheme_(fluxScheme::New(mesh))
{}


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
    if (oldIs_[stepi - 1] != -1)
    {
        rhoUOld_.set
        (
            oldIs_[stepi - 1],
            new volVectorField(rhoU_)
        );
        rhoEOld_.set
        (
            oldIs_[stepi - 1],
            new volScalarField(rhoE_)
        );
    }
    volVectorField rhoUOld(ai[stepi - 1]*rhoU_);
    volScalarField rhoEOld(ai[stepi - 1]*rhoE_);

    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = oldIs_[i];
        if (fi != -1 && ai[fi] != 0)
        {
            rhoUOld += ai[fi]*rhoUOld_[fi];
            rhoEOld += ai[fi]*rhoEOld_[fi];
        }
    }

    volVectorField deltaRhoU(fvc::div(rhoUPhi_));
    volScalarField deltaRhoE
    (
        fvc::div(rhoEPhi_)
      - ESource()
    );

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
    rhoU_ =
        cmptMultiply
        (
            rhoUOld
          - dT
           *(
                deltaRhoU
              - rho_*rho_.mesh().lookupObject<uniformDimensionedVectorField>("g")
            ),
            solutionDs
        );
    rhoE_ = rhoEOld - dT*deltaRhoE;
    if (radiation_->type() != "none")
    {
        calcAlphaAndRho();
        volScalarField T3(pow3(T()));

        volScalarField den(rho_ + f*dT*4.0*radiation_->Rp()*T3/Cv());

        volScalarField eNew
        (
            (rhoE_ - f*dT*radiation_->Rp()*T3*(T() - 4.0*e()/Cv()))/den
        );
        eNew.ref() += f*dT*radiation_->Ru()/den();
        rhoE_ = rho_*eNew;
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
         || UCoeff_.valid()
         || USource_.valid()
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
        if (UCoeff_.valid())
        {
            UEqn += fvm::Sp(UCoeff_(), U_);
        }
        if (USource_.valid())
        {
            UEqn += USource_();
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
            rhoE_ = rho_*(e() + 0.5*magSqr(U_)); // Includes change to total energy from viscos term in momentum equation

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

    UCoeff_.clear();
    USource_.clear();
}

void Foam::phaseCompressibleSystem::addUCoeff(const volScalarField& UCoeff)
{
    if (!UCoeff_.valid())
    {
        UCoeff_ = tmp<volScalarField>(new volScalarField("UCoeff", UCoeff));
    }
    else
    {
        UCoeff_.ref() += UCoeff;
    }
}


void Foam::phaseCompressibleSystem::addUSource(const volVectorField& USource)
{
    if (!USource_.valid())
    {
        USource_ = tmp<volVectorField>(new volVectorField("USource", USource));
    }
    else
    {
        USource_.ref() += USource;
    }
}


const Foam::fiveEqnCompressibleTurbulenceModel&
Foam::phaseCompressibleSystem::turbulence() const
{
    return turbulence_();
}


Foam::fiveEqnCompressibleTurbulenceModel&
Foam::phaseCompressibleSystem::turbulence()
{
    return turbulence_();
}


bool Foam::phaseCompressibleSystem::writeData(Ostream& os) const
{
    return os.good();
}

// ************************************************************************* //
