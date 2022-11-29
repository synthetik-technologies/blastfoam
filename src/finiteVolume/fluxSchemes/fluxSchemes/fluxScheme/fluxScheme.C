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

#include "fluxScheme.H"
#include "ReconstructionScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluxScheme, 0);
    defineRunTimeSelectionTable(fluxScheme, singlePhase);
    defineRunTimeSelectionTable(fluxScheme, multiphase);
    defineRunTimeSelectionTable(fluxScheme, interface);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxScheme::fluxScheme(const surfaceScalarField& phi)
:
    fluxSchemeBase(phi),
    dict_(mesh_.schemesDict())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxScheme::~fluxScheme()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluxScheme::clear()
{
    Uf_.clear();
}

void Foam::fluxScheme::createSavedFields()
{
    if (Uf_.valid())
    {
        return;
    }
    Uf_ = tmp<surfaceVectorField>
    (
        new surfaceVectorField
        (
            IOobject
            (
                fieldName("Uf"),
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedVector("0", dimVelocity, Zero)
        )
    );
}

Foam::tmp<Foam::surfaceVectorField> Foam::fluxScheme::Uf() const
{
    if (Uf_.valid())
    {
        return Uf_;
    }
    return surfaceVectorField::New
    (
       fieldName("Uf"),
        mesh_,
        dimensionedVector("0", dimVelocity, Zero)
    );
}


void Foam::fluxScheme::update
(
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& e,
    const volScalarField& p,
    const volScalarField& c,
    surfaceScalarField& phi,
    surfaceScalarField& rhoPhi,
    surfaceVectorField& rhoUPhi,
    surfaceScalarField& rhoEPhi
)
{
    createSavedFields();

    autoPtr<ReconstructionScheme<scalar>> rhoLimiter
    (
        ReconstructionScheme<scalar>::New(rho, "rho")
    );
    autoPtr<ReconstructionScheme<vector>> ULimiter
    (
        ReconstructionScheme<vector>::New(U, "U")
    );
    autoPtr<ReconstructionScheme<scalar>> eLimiter
    (
        ReconstructionScheme<scalar>::New(e, "e")
    );
    autoPtr<ReconstructionScheme<scalar>> pLimiter
    (
        ReconstructionScheme<scalar>::New(p, "p")
    );
    autoPtr<ReconstructionScheme<scalar>> cLimiter
    (
        ReconstructionScheme<scalar>::New(c, "speedOfSound")
    );

    tmp<surfaceScalarField> trhoOwn, trhoNei;
    rhoLimiter->interpolateOwnNei(trhoOwn, trhoNei);
    const surfaceScalarField& rhoOwn = trhoOwn();
    const surfaceScalarField& rhoNei = trhoNei();

    static bool cached = false;
    if (!cached)
    {
        cached = true;
        mesh_.addTemporaryObject(rhoOwn.name());
        mesh_.addTemporaryObject(rhoNei.name());
    }

    tmp<surfaceVectorField> tUOwn, tUNei;
    ULimiter->interpolateOwnNei(tUOwn, tUNei);
    const surfaceVectorField& UOwn = tUOwn();
    const surfaceVectorField& UNei = tUNei();

    tmp<surfaceScalarField> teOwn;
    tmp<surfaceScalarField> teNei;
    eLimiter->interpolateOwnNei(teOwn, teNei);
    const surfaceScalarField& eOwn = teOwn();
    const surfaceScalarField& eNei = teNei();

    tmp<surfaceScalarField> tpOwn;
    tmp<surfaceScalarField> tpNei;
    pLimiter->interpolateOwnNei(tpOwn, tpNei);
    const surfaceScalarField& pOwn = tpOwn();
    const surfaceScalarField& pNei = tpNei();

    tmp<surfaceScalarField> tcOwn;
    tmp<surfaceScalarField> tcNei;
    cLimiter->interpolateOwnNei(tcOwn, tcNei);
    const surfaceScalarField& cOwn = tcOwn();
    const surfaceScalarField& cNei = tcNei();

    preUpdate(p);
    forAll(UOwn, facei)
    {
        calculateFluxes
        (
            rhoOwn[facei], rhoNei[facei],
            UOwn[facei], UNei[facei],
            eOwn[facei], eNei[facei],
            pOwn[facei], pNei[facei],
            cOwn[facei], cNei[facei],
            mesh_.Sf()[facei],
            phi[facei],
            rhoPhi[facei],
            rhoUPhi[facei],
            rhoEPhi[facei],
            facei
        );
    }

    forAll(U.boundaryField(), patchi)
    {
        scalarField& pphi = phi.boundaryFieldRef()[patchi];
        scalarField& prhoPhi = rhoPhi.boundaryFieldRef()[patchi];
        vectorField& prhoUPhi = rhoUPhi.boundaryFieldRef()[patchi];
        scalarField& prhoEPhi = rhoEPhi.boundaryFieldRef()[patchi];
        forAll(U.boundaryField()[patchi], facei)
        {
            calculateFluxes
            (
                rhoOwn.boundaryField()[patchi][facei],
                rhoNei.boundaryField()[patchi][facei],
                UOwn.boundaryField()[patchi][facei],
                UNei.boundaryField()[patchi][facei],
                eOwn.boundaryField()[patchi][facei],
                eNei.boundaryField()[patchi][facei],
                pOwn.boundaryField()[patchi][facei],
                pNei.boundaryField()[patchi][facei],
                cOwn.boundaryField()[patchi][facei],
                cNei.boundaryField()[patchi][facei],
                mesh_.Sf().boundaryField()[patchi][facei],
                pphi[facei],
                prhoPhi[facei],
                prhoUPhi[facei],
                prhoEPhi[facei],
                facei, patchi
            );
        }
    }
    postUpdate();
}


Foam::tmp<Foam::surfaceScalarField> Foam::fluxScheme::energyFlux
(
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& e,
    const volScalarField& p
) const
{
    tmp<surfaceScalarField> trhoOwn, trhoNei;
    autoPtr<ReconstructionScheme<scalar>> rhoLimiter
    (
        ReconstructionScheme<scalar>::New(rho, "rho", false)
    );
    rhoLimiter->interpolateOwnNei(trhoOwn, trhoNei);
    const surfaceScalarField& rhoOwn = trhoOwn();
    const surfaceScalarField& rhoNei = trhoNei();

    // Interpolate fields
    autoPtr<ReconstructionScheme<vector>> ULimiter
    (
        ReconstructionScheme<vector>::New(U, "U", false)
    );
    autoPtr<ReconstructionScheme<scalar>> eLimiter
    (
        ReconstructionScheme<scalar>::New(e, "e")
    );
    autoPtr<ReconstructionScheme<scalar>> pLimiter
    (
        ReconstructionScheme<scalar>::New(p, "p", false)
    );

    tmp<surfaceVectorField> tUOwn, tUNei;
    ULimiter->interpolateOwnNei(tUOwn, tUNei);
    const surfaceVectorField& UOwn = tUOwn();
    const surfaceVectorField& UNei = tUNei();

    tmp<surfaceScalarField> teOwn, teNei;
    eLimiter->interpolateOwnNei(teOwn, teNei);
    const surfaceScalarField& eOwn = teOwn();
    const surfaceScalarField& eNei = teNei();

    tmp<surfaceScalarField> tpOwn, tpNei;
    pLimiter->interpolateOwnNei(tpOwn, tpNei);
    const surfaceScalarField& pOwn = tpOwn();
    const surfaceScalarField& pNei = tpNei();

    static bool cached = false;
    if (!cached)
    {
        cached = true;
        mesh_.addTemporaryObject(rhoOwn.name());
        mesh_.addTemporaryObject(rhoNei.name());
        mesh_.addTemporaryObject(UOwn.name());
        mesh_.addTemporaryObject(UNei.name());
        mesh_.addTemporaryObject(eOwn.name());
        mesh_.addTemporaryObject(eNei.name());
        mesh_.addTemporaryObject(pOwn.name());
        mesh_.addTemporaryObject(pNei.name());
    }

    tmp<surfaceScalarField> tmpPhi
    (
        new surfaceScalarField
        (
            IOobject
            (
                e.name() + "Phi",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar
            (
                "0",
                rho.dimensions()*e.dimensions()*dimVelocity*dimArea,
                0.0
            )
        )
    );
    surfaceScalarField& phi = tmpPhi.ref();

    forAll(eOwn, facei)
    {
        phi[facei] = energyFlux
        (
            rhoOwn[facei], rhoNei[facei],
            UOwn[facei], UNei[facei],
            eOwn[facei], eNei[facei],
            pOwn[facei], pNei[facei],
            mesh_.Sf()[facei],
            facei
        );
    }

    forAll(e.boundaryField(), patchi)
    {
        forAll(e.boundaryField()[patchi], facei)
        {
            phi.boundaryFieldRef()[patchi][facei] =
                energyFlux
                (
                    rhoOwn.boundaryField()[patchi][facei],
                    rhoNei.boundaryField()[patchi][facei],
                    UOwn.boundaryField()[patchi][facei],
                    UNei.boundaryField()[patchi][facei],
                    eOwn.boundaryField()[patchi][facei],
                    eNei.boundaryField()[patchi][facei],
                    pOwn.boundaryField()[patchi][facei],
                    pNei.boundaryField()[patchi][facei],
                    mesh_.Sf().boundaryField()[patchi][facei],
                    facei, patchi
                );
        }
    }
    return tmpPhi;
}

// ************************************************************************* //
