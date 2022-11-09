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

#include "shallowWaterSystem.H"
#include "fvm.H"
#include "hUInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static member functions * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(shallowWaterSystem, 0);

    template<>
    const char* Foam::NamedEnum<shallowWaterSystem::Friction, 3>::names[] =
    {
        "none",
        "Manning",
        "DarcyWeisbach"
    };

    const Foam::NamedEnum<shallowWaterSystem::Friction, 3>
        shallowWaterSystem::frictionTypes;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shallowWaterSystem::shallowWaterSystem
(
    const fvMesh& mesh
)
:
    timeIntegrationSystem(typeName, mesh),
    dict_
    (
        IOobject
        (
            "shallowWaterProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    g_
    (
        IOobject
        (
            "g",
            mesh.time().constant(),
            mesh
        ),
        dimensionedVector("g", dimAcceleration, dict_)
    ),
    rotating_(dict_.lookup<bool>("rotating")),
    omega_("omega", inv(dimTime), dict_.lookupOrDefault("omega", vector::zero)),
    F_("F", ((2.0*omega_ & g_)*g_/magSqr(g_))),

    h_
    (
        IOobject
        (
            "h",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    h0_
    (
        IOobject
        (
            "h0",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimLength, Zero),
        "zeroGradient"
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

    hU_
    (
        IOobject
        (
            "hU",
            mesh.time().timeName(),
            mesh
        ),
        h_*U_
    ),

    phi_
    (
        IOobject
        (
            "phi",
            mesh.time().timeName(),
            mesh
        ),
        fvc::flux(U_)
    ),

    hPhi_
    (
        IOobject
        (
            "hPhi",
            mesh.time().timeName(),
            mesh
        ),
        fvc::flux(hU_)
    ),

    hUPhi_
    (
        IOobject
        (
            "hUPhi",
            mesh.time().timeName(),
            mesh
        ),
        hPhi_*fvc::interpolate(U_)
    ),

    friction_(dict_.lookup<bool>("friction")),
    frictionType_
    (
        friction_
      ? frictionTypes.read(dict_.lookup("frictionType"))
      : none
    ),

    nPtr_(nullptr),
    fPtr_(nullptr),

    viscous_(dict_.lookupOrDefault<bool>("viscous", false)),
    turbulence_(viscous_ ? dict_.lookup<bool>("turbulence") : false),
    muh_("muh", dimensionSet(0, 2, -1, 0, 0, 0, 0), 0.0),
    muv_("muv", dimensionSet(0, 2, -1, 0, 0, 0, 0), 0.0),
    kl_("kl", dimensionSet(0, 1, -1, 0, 0, 0, 0), 0.0),
    kt_("kt", dimensionSet(0, -1, 0, 0, 0, 0, 0), 0.0),

    rain_(dict_.lookupOrDefault<bool>("rain", false)),
    rainfall_(nullptr),

    flux_(shallowFluxScheme::New(phi_, hPhi_, hUPhi_, g_)),
    hMin_("hMin", dimLength, dict_.lookupOrDefault<scalar>("hMin", 1e-6))
{

    h0_.correctBoundaryConditions();
    if (!mesh.time().restart())
    {
        IOobject hTotalHeader
        (
            "hTotal",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if (hTotalHeader.typeHeaderOk<volScalarField>())
        {
            volScalarField hTotal(hTotalHeader, mesh);
            h_ = hTotal - h0_;
            h_.correctBoundaryConditions();
            h_.write();
        }
    }

    switch (frictionType_)
    {
        case Manning:
        {
            IOobject nHeader
            (
                "n",
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            );

            if (!nHeader.typeHeaderOk<volScalarField>() && !dict_.found("n"))
            {
                FatalErrorInFunction
                    << "Friction is turned on, but " << string("n") << " was not" << nl
                    << "provided. " << endl
                    << abort(FatalError);
            }

            nPtr_.set
            (
                new volScalarField
                (
                    nHeader,
                    mesh,
                    dimensionedScalar
                    (
                        "n",
                        dimensionSet(0, -1.0/3.0, 1, 0, 0, 0, 0),
                        dict_.lookupOrDefault<scalar>("n", 1.0)
                    )
                )
            );
            nHeader.instance() = mesh.time().timeName();
            break;
        }
        case DarcyWeisbach:
        {
            IOobject fHeader
            (
                "f",
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            );

            if (!fHeader.typeHeaderOk<volScalarField>() && !dict_.found("f"))
            {
                FatalErrorInFunction
                    << "Friction is turned on, but " << string("f") << " was not" << nl
                    << "provided. " << endl
                    << abort(FatalError);
            }
            fPtr_.set
            (
                new volScalarField
                (
                    fHeader,
                    mesh,
                    dimensionedScalar
                    (
                        "f",
                        dimless,
                        dict_.lookupOrDefault<scalar>("f", 1.0)
                    )
                )
            );
            fHeader.instance() = mesh.time().timeName();
            break;
        }
        default:
        {
            break;
        }
    }

    if (viscous_)
    {
        if (dict_.found("mu"))
        {
            muh_.value() = dict_.lookup<scalar>("mu")/4.0;
        }
        else if (dict_.found("muh"))
        {
            muh_.read(dict_);
        }
        else
        {
            FatalIOErrorInFunction(dict_)
                << "Viscosity is turned on, but neither " << string("mu")
                << " or " << string("muh") << " was found" << endl
                << abort(FatalIOError);
        }

        if (turbulence_)
        {
            muv_.read(dict_);
            kl_.read(dict_);
            kt_.read(dict_);
        }
        else
        {
            kl_.readIfPresent(dict_);
        }
    }

    if (rain_)
    {
        rainfall_ = rainfallModel::New(mesh, dict_);
    }

    encode();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::shallowWaterSystem::~shallowWaterSystem()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::shallowWaterSystem::solve()
{
    volScalarField hDelta("hDelta", fvc::div(hPhi_));
    volVectorField hUDelta
    (
        "hUDelta",
        fvc::div(hUPhi_) + flux_->ghGradH0(g_, h_, h0_)
    );
    if (rotating_)
    {
        hUDelta += (F_ ^ hU_);
    }
    if (rain_)
    {
        hDelta -= rainfall_->R0();
    }

    //- Store changed in mass, momentum and energy
    this->storeAndBlendDelta(hDelta);
    this->storeAndBlendDelta(hUDelta);

    //- Store old values
    this->storeAndBlendOld(h_);
    this->storeAndBlendOld(hU_);

    dimensionedScalar dT = mesh().time().deltaT();
    h_ -= dT*hDelta;
    hU_ -= dT*hUDelta;
}


void Foam::shallowWaterSystem::postUpdate()
{
    this->decode();

    // Solve mass
    if (needSolve(h_.name()))
    {
        fvScalarMatrix hEqn
        (
            fvm::ddt(h_) - fvc::ddt(h_)
        ==
            models().source(h_)
        );

        constraints().constrain(hEqn);
        hEqn.solve();
        constraints().constrain(h_);
    }

    if (friction_ || viscous_)
    {
        volScalarField SfByU
        (
            volScalarField::New
            (
                "SfByU",
                mesh(),
                dimensionedScalar(inv(dimTime), 0.0)
            )
        );
        volVectorField UStar(U_);
        if (friction_)
        {
            tmp<volScalarField> K;
            if (frictionType_ == Manning)
            {
                K = sqr(nPtr_())/pow(max(h_, hMin_), 4.0/3.0);
            }
            else if (frictionType_ == DarcyWeisbach)
            {
                K = fPtr_()/((mag(g_)*8.0)*max(h_, hMin_));
            }
            SfByU += mag(g_)*K*mag(U_);
        }
        if (viscous_)
        {
            volScalarField hValid(pos(h_ - hMin_));
            volScalarField hs(max(h_, hMin_));
            if (kl_.value() > 0)
            {
                SfByU += kl_/(1.0 + kl_*h_/(3.0*muv_))*hValid/hs;
            }
            if (turbulence_)
            {
                SfByU += mag(U_)*kt_/sqr(1.0 + kl_*h_/(3.0*muv_));
            }
            UStar += mesh().time().deltaT()*(muh_*4.0)*fvc::laplacian(h_, U_)*hValid/hs;
        }

        U_ = UStar/(1.0 + mesh().time().deltaT()*SfByU);
        U_.correctBoundaryConditions();
        hU_ = h_*U_;
    }

    if (needSolve(U_.name()))
    {
        fvVectorMatrix UEqn
        (
            fvm::ddt(h_, U_) - fvc::ddt(hU_)
          + fvm::ddt(hMin_, U_) - fvc::ddt(hMin_, U_)
        ==
            models().source(h_, U_)
        );

        constraints().constrain(UEqn);
        UEqn.solve();
        constraints().constrain(U_);

        hU_ = h_*U_;
    }
}


void Foam::shallowWaterSystem::update()
{
    decode();
    flux_->update
    (
        h_,
        h0_,
        U_,
        hU_
    );
}


void Foam::shallowWaterSystem::decode()
{
    h_.max(0);
    h_.correctBoundaryConditions();

    U_.ref() = hU_()/max(h_(), hMin_);
    U_.correctBoundaryConditions();
    hU_ = h_*U_;
}


void Foam::shallowWaterSystem::encode()
{
    hU_ = h_*U_;
}


Foam::scalar Foam::shallowWaterSystem::CoNum() const
{
    const fvMesh& mesh = this->mesh();
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;
    scalar waveCoNum = 0.0;

    if (mesh.nInternalFaces())
    {
        const scalarField& V = mesh.V();
        const scalar& deltaT = mesh.time().deltaTValue();

        surfaceScalarField ws(sqrt(fvc::interpolate(h_)*mag(g_)));

        scalarField sumPhi
        (
            fvc::surfaceSum(mag(phi_) + ws*mesh.magSf()*0.5)().primitiveField()
        );

        CoNum = 0.5*gMax(sumPhi/V)*deltaT;

        meanCoNum = 0.5*(gSum(sumPhi)/gSum(V))*deltaT;

        // Gravity wave Courant number
        waveCoNum = 0.25*gMax
        (
            fvc::surfaceSum(ws*mesh.magSf())().primitiveField()/V
        )*deltaT;
    }

    Info<< "Courant number mean: " << meanCoNum
        << " max: " << CoNum << endl;

    Info<< "Gravity wave Courant number max: " << waveCoNum
        << endl;

    return CoNum;
}

// ************************************************************************* //
