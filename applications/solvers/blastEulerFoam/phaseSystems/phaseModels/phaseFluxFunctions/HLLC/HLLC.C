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

#include "HLLC.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluxSchemes
{
    defineTypeNameAndDebug(HLLC, 0);
    addToRunTimeSelectionTable(fluxScheme, HLLC, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxSchemes::HLLC::HLLC
(
    const fvMesh& mesh,
    const word& name
)
:
    fluxScheme(mesh, name)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxSchemes::HLLC::~HLLC()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluxSchemes::HLLC::clear()
{
    fluxScheme::clear();
    SOwn_.clear();
    SNei_.clear();
    SStar_.clear();
    UvOwn_.clear();
    UvNei_.clear();
}

void Foam::fluxSchemes::HLLC::createSavedFields()
{
    fluxScheme::createSavedFields();
    if (SOwn_.valid())
    {
        return;
    }
    SOwn_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                IOobject::groupName("HLLC::SOwn", this->group()),
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimVelocity, 0.0)
        )
    );
    SNei_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                IOobject::groupName("HLLC::SNei", this->group()),
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimVelocity, 0.0)
        )
    );
    SStar_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                IOobject::groupName("HLLC::SStar", this->group()),
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimVelocity, 0.0)
        )
    );
    UvOwn_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                IOobject::groupName("HLLC::UvOwn", this->group()),
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimVelocity, 0.0)
        )
    );
    UvNei_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                IOobject::groupName("HLLC::UvNei", this->group()),
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimVelocity, 0.0)
        )
    );
}


void Foam::fluxSchemes::HLLC::calculateFluxes
(
    const scalar& alphaOwn, const scalar& alphaNei,
    const scalar& rhoO, const scalar& rhoN,
    const vector& UOwn, const vector& UNei,
    const scalar& eOwn, const scalar& eNei,
    const scalar& pOwn, const scalar& pNei,
    const scalar& cOwn, const scalar& cNei,
    const vector& Sf,
    scalar& phi,
    scalar& alphaPhi,
    scalar& alphaRhoPhi,
    vector& alphaRhoUPhi,
    scalar& alphaRhoEPhi,
    const label facei, const label patchi
)
{
    scalar magSf = mag(Sf);
    vector normal = Sf/magSf;

    scalar EOwn = eOwn + 0.5*magSqr(UOwn);
    scalar ENei = eNei + 0.5*magSqr(UNei);

    const scalar vMesh(meshPhi(facei, patchi)/magSf);
    scalar UvOwn((UOwn & normal) - vMesh);
    scalar UvNei((UNei & normal) - vMesh);

    scalar rhoOwn = max(rhoO, 1e-10);
    scalar rhoNei = max(rhoN, 1e-10);

    scalar wOwn(sqrt(rhoOwn)/(sqrt(rhoOwn) + sqrt(rhoNei)));
    scalar wNei(1.0 - wOwn);

    scalar cTilde(cOwn*wOwn + cNei*wNei);
    vector UTilde(UOwn*wOwn + UNei*wNei);
    scalar UvTilde(UTilde & normal);

    scalar SOwn(min(UvOwn - cOwn, UvTilde - cTilde));
    scalar SNei(max(UvNei + cNei, UvTilde + cTilde));

    scalar SStar
    (
        (
            pNei - pOwn
          + rhoOwn*UvOwn*(SOwn - UvOwn)
          - rhoNei*UvNei*(SNei - UvNei)
        )
       /stabilise(rhoOwn*(SOwn - UvOwn) - rhoNei*(SNei - UvNei), small)
    );

    scalar pStarOwn(pOwn + rhoOwn*(SOwn - UvOwn)*(SStar - UvOwn));
    scalar pStarNei(pNei + rhoNei*(SNei - UvNei)*(SStar - UvNei));
    scalar pStar(0.5*(pStarOwn + pStarNei));

    this->save(facei, patchi, SOwn, SOwn_);
    this->save(facei, patchi, SNei, SNei_);
    this->save(facei, patchi, SStar, SStar_);
    this->save(facei, patchi, UvOwn, UvOwn_);
    this->save(facei, patchi, UvNei, UvNei_);

    // Owner values
    scalar alpha;
    scalar rho;
    vector U;
    scalar E;
    scalar p;

    if (SOwn > 0)
    {
        alpha = alphaOwn;
        rho = rhoOwn;
        phi = UvOwn;
        U = UOwn;
        E = EOwn;
        p = pOwn;
    }
    else if (SStar > 0)
    {
        alpha = alphaOwn;
        rho = rhoOwn*(SOwn - UvOwn)/(SOwn - SStar);
        phi = SStar;
        U = (UOwn - UvOwn*normal) + SStar*normal;
        E = EOwn + (pStar*SStar - pOwn*UvOwn)/(rhoOwn*(SOwn - UvOwn));
        p = pStar;
    }
    else if (SNei > 0)
    {
        alpha = alphaNei;
        rho = rhoNei*(SNei - UvNei)/(SNei - SStar);
        phi = SStar;
        U = (UNei - UvNei*normal) + SStar*normal;
        E = ENei + (pStar*SStar - pNei*UvNei)/(rhoNei*(SNei - UvNei));
        p = pStar;
    }
    else
    {
        alpha = alphaNei;
        rho = rhoNei;
        phi = UvNei;
        U = UNei;
        E = ENei;
        p = pNei;
    }

    this->save(facei, patchi, alpha, alphaf_);
    this->save(facei, patchi, U, Uf_);
    this->save(facei, patchi, p, pf_);

    phi *= magSf;
    alphaPhi = alpha*phi;
    alphaRhoPhi = rho*alphaPhi;
    alphaRhoUPhi = alphaRhoPhi*U + alpha*p*Sf;
    alphaRhoEPhi = alphaPhi*(rho*E + p);
    alphaRhoEPhi += vMesh*magSf*p*alpha;
}


void Foam::fluxSchemes::HLLC::calculateFluxes
(
    const scalar& alphaOwn, const scalar& alphaNei,
    const scalar& rhoO, const scalar& rhoN,
    const scalarList& alphasOwn, const scalarList& alphasNei,
    const scalarList& rhosOwn, const scalarList& rhosNei,
    const vector& UOwn, const vector& UNei,
    const scalar& eOwn, const scalar& eNei,
    const scalar& pOwn, const scalar& pNei,
    const scalar& cOwn, const scalar& cNei,
    const vector& Sf,
    scalar& phi,
    scalarList& alphaPhis,
    scalarList& alphaRhoPhis,
    vector& alphaRhoUPhi,
    scalar& alphaRhoEPhi,
    const label facei, const label patchi
)
{
    scalar magSf = mag(Sf);
    vector normal = Sf/magSf;

    scalar EOwn = eOwn + 0.5*magSqr(UOwn);
    scalar ENei = eNei + 0.5*magSqr(UNei);

    const scalar vMesh(meshPhi(facei, patchi)/magSf);
    scalar UvOwn((UOwn & normal) - vMesh);
    scalar UvNei((UNei & normal) - vMesh);

    scalar rhoOwn = max(rhoO, 1e-10);
    scalar rhoNei = max(rhoN, 1e-10);

    scalar wOwn(sqrt(rhoOwn)/(sqrt(rhoOwn) + sqrt(rhoNei)));
    scalar wNei(1.0 - wOwn);

    scalar cTilde(cOwn*wOwn + cNei*wNei);
    vector UTilde(UOwn*wOwn + UNei*wNei);
    scalar UvTilde(UTilde & normal);

    scalar SOwn(min(UvOwn - cOwn, UvTilde - cTilde));
    scalar SNei(max(UvNei + cNei, UvTilde + cTilde));

    scalar SStar
    (
        (
            pNei - pOwn
          + rhoOwn*UvOwn*(SOwn - UvOwn)
          - rhoNei*UvNei*(SNei - UvNei)
        )
       /stabilise(rhoOwn*(SOwn - UvOwn) - rhoNei*(SNei - UvNei), small)
    );

    scalar pStarOwn(pOwn + rhoOwn*(SOwn - UvOwn)*(SStar - UvOwn));
    scalar pStarNei(pNei + rhoNei*(SNei - UvNei)*(SStar - UvNei));
    scalar pStar(0.5*(pStarOwn + pStarNei));

    this->save(facei, patchi, SOwn, SOwn_);
    this->save(facei, patchi, SNei, SNei_);
    this->save(facei, patchi, SStar, SStar_);
    this->save(facei, patchi, UvOwn, UvOwn_);
    this->save(facei, patchi, UvNei, UvNei_);

    // Owner values
    scalar alpha;
    scalar rho;
    scalarList alphas(alphasOwn.size());
    scalarList rhos(rhosOwn.size());
    vector U;
    scalar E;
    scalar p;

    if (SOwn > 0)
    {
        alpha = alphaOwn;
        rho = rhoOwn;
        phi = UvOwn;
        U = UOwn;
        E = EOwn;
        p = pOwn;

        forAll(alphas, phasei)
        {
            alphas[phasei] = alphasOwn[phasei];
            rhos[phasei] = rhosOwn[phasei];
        }
    }
    else if (SStar > 0)
    {
        scalar f = (SOwn - UvOwn)/(SOwn - SStar);
        alpha = alphaOwn;
        rho = rhoOwn*f;
        phi = SStar;
        U = (UOwn - UvOwn*normal) + SStar*normal;
        E = EOwn + (pStar*SStar - pOwn*UvOwn)/(rhoOwn*(SOwn - UvOwn));
        p = pStar;

        forAll(alphas, phasei)
        {
            alphas[phasei] = alphasOwn[phasei];
            rhos[phasei] = rhosOwn[phasei]*f;
        }
    }
    else if (SNei > 0)
    {
        scalar f = (SNei - UvNei)/(SNei - SStar);
        alpha = alphaNei;
        rho = rhoNei*f;
        phi = SStar;
        U = (UNei - UvNei*normal) + SStar*normal;
        E = ENei + (pStar*SStar - pNei*UvNei)/(rhoNei*(SNei - UvNei));
        p = pStar;

        forAll(alphas, phasei)
        {
            alphas[phasei] = alphasNei[phasei];
            rhos[phasei] = rhosNei[phasei]*f;
        }
    }
    else
    {
        alpha = alphaNei;
        rho = rhoNei;
        phi = UvNei;
        U = UNei;
        E = ENei;
        p = pNei;

        forAll(alphas, phasei)
        {
            alphas[phasei] = alphasNei[phasei];
            rhos[phasei] = rhosNei[phasei];
        }
    }

    this->save(facei, patchi, alpha, alphaf_);
    this->save(facei, patchi, U, Uf_);
    this->save(facei, patchi, p, pf_);

    phi *= magSf;
    alphaRhoUPhi = alpha*(rho*U*phi + p*Sf);
    alphaRhoEPhi = alpha*phi*(rho*E + p);
    alphaRhoEPhi += vMesh*magSf*p*alpha;

    forAll(alphaPhis, phasei)
    {
        alphaPhis[phasei] = alphas[phasei]*phi;
        alphaRhoPhis[phasei] = alphaPhis[phasei]*rhos[phasei];
    }
}


Foam::scalar Foam::fluxSchemes::HLLC::interpolate
(
    const scalar& fOwn, const scalar& fNei,
    const bool rho,
    const label facei, const label patchi
) const
{
    scalar SOwn = getValue(facei, patchi, SOwn_());
    scalar SStar = getValue(facei, patchi, SStar_());

    if (rho)
    {
        scalar SNei = getValue(facei, patchi, SNei_());
        scalar UvOwn = getValue(facei, patchi, UvOwn_());
        scalar UvNei = getValue(facei, patchi, UvNei_());
        if (SOwn > 0)
        {
            return fOwn;
        }
        else if (SStar > 0)
        {
            return fOwn*(SOwn - UvOwn)/(SOwn - SStar);
        }
        else if (SNei > 0)
        {
            return fNei*(SNei - UvNei)/(SNei - SStar);
        }
        return fNei;
    }
    if (SOwn > 0 || SStar > 0)
    {
        return fOwn;
    }
    else
    {
        return fNei;
    }
}

// ************************************************************************* //
