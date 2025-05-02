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

#include "phaseFluxScheme.H"
#include "ReconstructionScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseFluxScheme, 0);
    defineRunTimeSelectionTable(phaseFluxScheme, dictionary);
    defineRunTimeSelectionTable(phaseFluxScheme, solid);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseFluxScheme::phaseFluxScheme(const surfaceScalarField& phi)
:
    fluxSchemeBase(phi),
    phaseName_(phi.group()),
    dict_
    (
        phi.mesh().schemes().dict().subDict
        (
            "fluxSchemes"
        ).subDict(phaseName_)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseFluxScheme::~phaseFluxScheme()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::phaseFluxScheme::clear()
{
    Uf_.clear();
    pf_.clear();
    alphaf_.clear();
}

void Foam::phaseFluxScheme::createSavedFields()
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
                mesh_.time().name(),
                mesh_
            ),
            mesh_,
            dimensionedVector("0", dimVelocity, Zero)
        )
    );
    pf_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                fieldName("pf"),
                mesh_.time().name(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimPressure, Zero)
        )
    );
    alphaf_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                fieldName("alphaf"),
                mesh_.time().name(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimless, Zero)
        )
    );
}

Foam::tmp<Foam::surfaceVectorField> Foam::phaseFluxScheme::Uf() const
{
    if (Uf_.valid())
    {
        return Uf_();
    }
    FatalErrorInFunction
        << fieldName("Uf")
        << " has not been set." << nl
        << abort(FatalError);

    return Uf_();
}


Foam::tmp<Foam::surfaceScalarField> Foam::phaseFluxScheme::pf() const
{
    if (pf_.valid())
    {
        return pf_();
    }
    FatalErrorInFunction
        << fieldName("pf")
        << " has not been set." << nl
        << abort(FatalError);

    return pf_;
}


Foam::tmp<Foam::surfaceScalarField> Foam::phaseFluxScheme::alphaf() const
{
    if (alphaf_.valid())
    {
        return alphaf_();
    }
    FatalErrorInFunction
        << fieldName("alphaf")
        << " has not been set." << nl
        << abort(FatalError);

    return alphaf_;
}


void Foam::phaseFluxScheme::update
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& e,
    const volScalarField& p,
    const volScalarField& c,
    surfaceScalarField& phi,
    surfaceScalarField& alphaPhi,
    surfaceScalarField& alphaRhoPhi,
    surfaceVectorField& alphaRhoUPhi,
    surfaceScalarField& alphaRhoEPhi
)
{
    createSavedFields();

    autoPtr<ReconstructionScheme<scalar>> alphaLimiter
    (
        ReconstructionScheme<scalar>::New(alpha, "alpha", phaseName_, true)
    );
    autoPtr<ReconstructionScheme<scalar>> rhoLimiter
    (
        ReconstructionScheme<scalar>::New(rho, "rho", phaseName_, true)
    );
    autoPtr<ReconstructionScheme<vector>> ULimiter
    (
        ReconstructionScheme<vector>::New(U, "U", phaseName_, true)
    );
    autoPtr<ReconstructionScheme<scalar>> eLimiter
    (
        ReconstructionScheme<scalar>::New(e, "e", phaseName_, true)
    );
    autoPtr<ReconstructionScheme<scalar>> pLimiter
    (
        ReconstructionScheme<scalar>::New(p, "p", phaseName_, true)
    );
    autoPtr<ReconstructionScheme<scalar>> cLimiter
    (
        ReconstructionScheme<scalar>::New(c, "speedOfSound", phaseName_, true)
    );

    tmp<surfaceScalarField> talphaOwn;
    tmp<surfaceScalarField> talphaNei;
    alphaLimiter->interpolateOwnNei(talphaOwn, talphaNei);
    const surfaceScalarField& alphaOwn = talphaOwn();
    const surfaceScalarField& alphaNei = talphaNei();

    tmp<surfaceScalarField> trhoOwn;
    tmp<surfaceScalarField> trhoNei;
    rhoLimiter->interpolateOwnNei(trhoOwn, trhoNei);
    const surfaceScalarField& rhoOwn = trhoOwn();
    const surfaceScalarField& rhoNei = trhoNei();

    if
    (
        mesh_.cacheTemporaryObject
        (
            reconstruction::ownName(IOobject::groupName("alphaRho", phaseName_))
        )
    )
    {
        surfaceScalarField alphaRhoNei
        (
            reconstruction::ownName(IOobject::groupName("alphaRho", phaseName_)),
            alphaOwn*rhoOwn
        );
        mesh_.cacheTemporaryObject(alphaRhoNei);
    }
    if
    (
        mesh_.cacheTemporaryObject
        (
            reconstruction::neiName(IOobject::groupName("alphaRho", phaseName_))
        )
    )
    {
        surfaceScalarField alphaRhoNei
        (
            reconstruction::neiName(IOobject::groupName("alphaRho", phaseName_)),
            alphaNei*rhoNei
        );
        mesh_.cacheTemporaryObject(alphaRhoNei);
    }

    tmp<surfaceVectorField> tUOwn;
    tmp<surfaceVectorField> tUNei;
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
            alphaOwn[facei], alphaNei[facei],
            rhoOwn[facei], rhoNei[facei],
            UOwn[facei], UNei[facei],
            eOwn[facei], eNei[facei],
            pOwn[facei], pNei[facei],
            cOwn[facei], cNei[facei],
            mesh_.Sf()[facei],
            phi[facei],
            alphaPhi[facei],
            alphaRhoPhi[facei],
            alphaRhoUPhi[facei],
            alphaRhoEPhi[facei],
            facei
        );
    }

    forAll(U.boundaryField(), patchi)
    {
        forAll(U.boundaryField()[patchi], facei)
        {
            calculateFluxes
            (
                alphaOwn.boundaryField()[patchi][facei],
                alphaNei.boundaryField()[patchi][facei],
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
                phi.boundaryFieldRef()[patchi][facei],
                alphaPhi.boundaryFieldRef()[patchi][facei],
                alphaRhoPhi.boundaryFieldRef()[patchi][facei],
                alphaRhoUPhi.boundaryFieldRef()[patchi][facei],
                alphaRhoEPhi.boundaryFieldRef()[patchi][facei],
                facei, patchi
            );
        }
    }
    postUpdate();
}

void Foam::phaseFluxScheme::update
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& e,
    const volScalarField& p,
    const volScalarField& c,
    surfaceScalarField& phi,
    surfaceScalarField& alphaRhoPhi,
    surfaceVectorField& alphaRhoUPhi,
    surfaceScalarField& alphaRhoEPhi
)
{
    createSavedFields();

    autoPtr<ReconstructionScheme<scalar>> alphaLimiter
    (
        ReconstructionScheme<scalar>::New(alpha, "alpha", phaseName_, true)
    );
    autoPtr<ReconstructionScheme<scalar>> rhoLimiter
    (
        ReconstructionScheme<scalar>::New(rho, "rho", phaseName_, true)
    );
    autoPtr<ReconstructionScheme<vector>> ULimiter
    (
        ReconstructionScheme<vector>::New(U, "U", phaseName_, true)
    );
    autoPtr<ReconstructionScheme<scalar>> eLimiter
    (
        ReconstructionScheme<scalar>::New(e, "e", phaseName_, true)
    );
    autoPtr<ReconstructionScheme<scalar>> pLimiter
    (
        ReconstructionScheme<scalar>::New(p, "p", phaseName_, true)
    );
    autoPtr<ReconstructionScheme<scalar>> cLimiter
    (
        ReconstructionScheme<scalar>::New(c, "speedOfSound", phaseName_, true)
    );

    tmp<surfaceScalarField> talphaOwn;
    tmp<surfaceScalarField> talphaNei;
    alphaLimiter->interpolateOwnNei(talphaOwn, talphaNei);
    const surfaceScalarField& alphaOwn = talphaOwn();
    const surfaceScalarField& alphaNei = talphaNei();

    tmp<surfaceScalarField> trhoOwn;
    tmp<surfaceScalarField> trhoNei;
    rhoLimiter->interpolateOwnNei(trhoOwn, trhoNei);
    const surfaceScalarField& rhoOwn = trhoOwn();
    const surfaceScalarField& rhoNei = trhoNei();

    if
    (
        mesh_.cacheTemporaryObject
        (
            reconstruction::ownName(IOobject::groupName("alphaRho", phaseName_))
        )
    )
    {
        surfaceScalarField alphaRhoNei
        (
            reconstruction::ownName(IOobject::groupName("alphaRho", phaseName_)),
            alphaOwn*rhoOwn
        );
        mesh_.cacheTemporaryObject(alphaRhoNei);
    }
    if
    (
        mesh_.cacheTemporaryObject
        (
            reconstruction::neiName(IOobject::groupName("alphaRho", phaseName_))
        )
    )
    {
        surfaceScalarField alphaRhoNei
        (
            reconstruction::neiName(IOobject::groupName("alphaRho", phaseName_)),
            alphaNei*rhoNei
        );
        mesh_.cacheTemporaryObject(alphaRhoNei);
    }

    tmp<surfaceVectorField> tUOwn;
    tmp<surfaceVectorField> tUNei;
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

    scalar alphaPhi;
    preUpdate(p);
    forAll(UOwn, facei)
    {
        calculateFluxes
        (
            alphaOwn[facei], alphaNei[facei],
            rhoOwn[facei], rhoNei[facei],
            UOwn[facei], UNei[facei],
            eOwn[facei], eNei[facei],
            pOwn[facei], pNei[facei],
            cOwn[facei], cNei[facei],
            mesh_.Sf()[facei],
            phi[facei],
            alphaPhi,
            alphaRhoPhi[facei],
            alphaRhoUPhi[facei],
            alphaRhoEPhi[facei],
            facei
        );
    }

    forAll(U.boundaryField(), patchi)
    {
        forAll(U.boundaryField()[patchi], facei)
        {
            calculateFluxes
            (
                alphaOwn.boundaryField()[patchi][facei],
                alphaNei.boundaryField()[patchi][facei],
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
                phi.boundaryFieldRef()[patchi][facei],
                alphaPhi,
                alphaRhoPhi.boundaryFieldRef()[patchi][facei],
                alphaRhoUPhi.boundaryFieldRef()[patchi][facei],
                alphaRhoEPhi.boundaryFieldRef()[patchi][facei],
                facei, patchi
            );
        }
    }
    postUpdate();
}

// ************************************************************************* //
