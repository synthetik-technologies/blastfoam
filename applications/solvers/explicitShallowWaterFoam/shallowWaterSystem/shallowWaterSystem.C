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
#include "triSurfaceMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static member functions * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(shallowWaterSystem, 0);
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
    nPtr_(nullptr),

    flux_(shallowFluxScheme::New(phi_, hPhi_, hUPhi_, g_)),
    hMin_("hMin", dimLength, dict_.lookupOrDefault<scalar>("hMin", 1e-6))
{

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

    if (friction_)
    {
        IOobject nHeader
        (
            "n",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        );

        if (!nHeader.typeHeaderOk<volScalarField>())
        {
            nHeader.instance() = mesh.time().timeName();
        }

        if (nHeader.typeHeaderOk<volScalarField>())
        {}
        else if (dict_.found("n"))
        {
            nHeader.readOpt() = IOobject::NO_READ;
        }
        else
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
        fvc::div(hUPhi_) + mag(g_)*h_*fvc::grad(h0_)
    );
    if (rotating_)
    {
        hUDelta += (F_ ^ hU_);
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

    if (friction_)
    {
        volScalarField K
        (
            volScalarField::New
            (
                "K",
                sqr(nPtr_())*mag(U_)/pow(max(h_, hMin_), 4.0/3.0)
            )
        );
        U_ = U_/(1.0 + mag(g_)*mesh().time().deltaT()*K);
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
        U_
    );
}


void Foam::shallowWaterSystem::decode()
{
    h_.max(0);
    h_.correctBoundaryConditions();

    U_.ref() = hU_()/max(h_(), hMin_);
    U_.correctBoundaryConditions();
    hU_ = h_*U_;
    // hU_.boundaryFieldRef() = U_.boundaryField()*h_.boundaryField();
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
