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
    const fvMesh& mesh
)
:
    fluxScheme(mesh)
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
    pStarOwn_.clear();
    pStarNei_.clear();
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
                "HLLC::SOwn",
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
                "HLLC::SNei",
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
                "HLLC::SStar",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimVelocity, 0.0)
        )
    );
    pStarOwn_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "HLLC::pStarOwn",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimPressure, 0.0)
        )
    );
    pStarNei_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "HLLC::pStarNei",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimPressure, 0.0)
        )
    );
    UvOwn_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "HLLC::UvOwn",
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
                "HLLC::UvNei",
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
    const scalar& rhoOwn, const scalar& rhoNei,
    const vector& UOwn, const vector& UNei,
    const scalar& eOwn, const scalar& eNei,
    const scalar& pOwn, const scalar& pNei,
    const scalar& cOwn, const scalar& cNei,
    const vector& Sf,
    scalar& phi,
    scalar& rhoPhi,
    vector& rhoUPhi,
    scalar& rhoEPhi,
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

    scalar wOwn(sqrt(rhoOwn)/(sqrt(rhoOwn) + sqrt(rhoNei)));
    scalar wNei(1.0 - wOwn);

    scalar cTilde(cOwn*wOwn + cNei*wNei);
    scalar UvTilde(UvOwn*wOwn + UvNei*wNei);

    scalar SOwn(min(UvOwn - cOwn, UvTilde - cTilde));
    scalar SNei(max(UvNei + cNei, UvTilde + cTilde));

    scalar SStar
    (
        (
            pNei - pOwn
          + rhoOwn*UvOwn*(SOwn - UvOwn)
          - rhoNei*UvNei*(SNei - UvNei)
        )
       /(rhoOwn*(SOwn - UvOwn) - rhoNei*(SNei - UvNei))
    );

    scalar pStarOwn(pOwn + rhoOwn*(SOwn - UvOwn)*(SStar - UvOwn));
    scalar pStarNei(pNei + rhoNei*(SNei - UvNei)*(SStar - UvNei));

    this->save(facei, patchi, SOwn, SOwn_);
    this->save(facei, patchi, SNei, SNei_);
    this->save(facei, patchi, SStar, SStar_);
    this->save(facei, patchi, pStarOwn, pStarOwn_);
    this->save(facei, patchi, pStarNei, pStarNei_);
    this->save(facei, patchi, UvOwn, UvOwn_);
    this->save(facei, patchi, UvNei, UvNei_);

    // Owner values
    const vector rhoUOwn = rhoOwn*UOwn;
    const scalar rhoEOwn = rhoOwn*EOwn;

    const vector rhoUPhiOwn = rhoUOwn*UvOwn + pOwn*normal;
    const scalar rhoEPhiOwn = (rhoEOwn + pOwn)*UvOwn;

    // Neighbour values
    const vector rhoUNei = rhoNei*UNei;
    const scalar rhoENei = rhoNei*ENei;

    const vector rhoUPhiNei = rhoUNei*UvNei + pNei*normal;
    const scalar rhoEPhiNei = (rhoENei + pNei)*UvNei;

    scalar p;
    if (SOwn > 0)
    {
        this->save(facei, patchi, UOwn, Uf_);
        phi = UvOwn;
        rhoPhi = rhoOwn*UvOwn;
        p = pOwn;
        rhoUPhi = rhoUPhiOwn;
        rhoEPhi = rhoEPhiOwn;
    }
    else if (SStar > 0)
    {
        const scalar dS = SOwn - SStar;

        this->save
        (
            facei,
            patchi,
            (SOwn*rhoUOwn - rhoUPhiOwn + pStarOwn*normal)
           /(rhoOwn*(SOwn - UvOwn)),
            Uf_
        );
        phi = SStar;
        rhoPhi = phi*rhoOwn*(SOwn - UvOwn)/dS;
        p = 0.5*(pStarNei + pStarOwn);
        rhoUPhi =
            (SStar*(SOwn*rhoUOwn - rhoUPhiOwn) + SOwn*pStarOwn*normal)/dS;
        rhoEPhi = SStar*(SOwn*rhoEOwn - rhoEPhiOwn + SOwn*pStarOwn)/dS;
    }
    else if (SNei > 0)
    {
        const scalar dS = SNei - SStar;

        this->save
        (
            facei,
            patchi,
            (SNei*rhoUNei - rhoUPhiNei + pStarNei*normal)
           /(rhoNei*(SNei - UvNei)),
            Uf_
        );
        phi = SStar;
        rhoPhi = phi*rhoNei*(SNei - UvNei)/dS;
        p = 0.5*(pStarNei + pStarOwn);
        rhoUPhi =
            (SStar*(SNei*rhoUNei - rhoUPhiNei) + SNei*pStarNei*normal)/dS;
        rhoEPhi = SStar*(SNei*rhoENei - rhoEPhiNei + SNei*pStarNei)/dS;
    }
    else
    {
        this->save(facei, patchi, UNei, Uf_);
        phi = UvNei;
        rhoPhi = rhoNei*UvNei;
        p = pNei;
        rhoUPhi = rhoUPhiNei;
        rhoEPhi = rhoEPhiNei;
    }
    phi *= magSf;
    rhoPhi *= magSf;
    rhoUPhi *= magSf;
    rhoEPhi *= magSf;
    rhoEPhi += vMesh*magSf*p;
}


void Foam::fluxSchemes::HLLC::calculateFluxes
(
    const scalarList& alphasOwn, const scalarList& alphasNei,
    const scalarList& rhosOwn, const scalarList& rhosNei,
    const scalar& rhoOwn, const scalar& rhoNei,
    const vector& UOwn, const vector& UNei,
    const scalar& eOwn, const scalar& eNei,
    const scalar& pOwn, const scalar& pNei,
    const scalar& cOwn, const scalar& cNei,
    const vector& Sf,
    scalar& phi,
    scalarList& alphaPhis,
    scalarList& alphaRhoPhis,
    vector& rhoUPhi,
    scalar& rhoEPhi,
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

    scalar wOwn(sqrt(rhoOwn)/(sqrt(rhoOwn) + sqrt(rhoNei)));
    scalar wNei(1.0 - wOwn);

    scalar cTilde(cOwn*wOwn + cNei*wNei);
    scalar UvTilde(UvOwn*wOwn + UvNei*wNei);

    scalar SOwn(min(UvOwn - cOwn, UvTilde - cTilde));
    scalar SNei(max(UvNei + cNei, UvTilde + cTilde));

    scalar SStar
    (
        (
            pNei - pOwn
          + rhoOwn*UvOwn*(SOwn - UvOwn)
          - rhoNei*UvNei*(SNei - UvNei)
        )
       /(rhoOwn*(SOwn - UvOwn) - rhoNei*(SNei - UvNei))
    );

    scalar pStarOwn(pOwn + rhoOwn*(SOwn - UvOwn)*(SStar - UvOwn));
    scalar pStarNei(pNei + rhoNei*(SNei - UvNei)*(SStar - UvNei));

    this->save(facei, patchi, SOwn, SOwn_);
    this->save(facei, patchi, SNei, SNei_);
    this->save(facei, patchi, SStar, SStar_);
    this->save(facei, patchi, pStarOwn, pStarOwn_);
    this->save(facei, patchi, pStarNei, pStarNei_);
    this->save(facei, patchi, UvOwn, UvOwn_);
    this->save(facei, patchi, UvNei, UvNei_);

    // Owner values
    const vector rhoUOwn = rhoOwn*UOwn;
    const scalar rhoEOwn = rhoOwn*EOwn;

    const vector rhoUPhiOwn = rhoUOwn*UvOwn + pOwn*normal;
    const scalar rhoEPhiOwn = (rhoEOwn + pOwn)*UvOwn;

    // Neighbour values
    const vector rhoUNei = rhoNei*UNei;
    const scalar rhoENei = rhoNei*ENei;

    const vector rhoUPhiNei = rhoUNei*UvNei + pNei*normal;
    const scalar rhoEPhiNei = (rhoENei + pNei)*UvNei;

    scalar p;
    if (SOwn > 0)
    {
        this->save(facei, patchi, UOwn, Uf_);
        phi = UvOwn;
        rhoUPhi = rhoUPhiOwn;
        p = pOwn;
        rhoEPhi = rhoEPhiOwn;

        forAll(alphaPhis, phasei)
        {
            alphaPhis[phasei] = alphasOwn[phasei]*UvOwn;
            alphaRhoPhis[phasei] = alphasOwn[phasei]*rhosOwn[phasei]*UvOwn;
        }
    }
    else if (SStar > 0)
    {
        const scalar dS = SOwn - SStar;

        this->save
        (
            facei,
            patchi,
            (SOwn*rhoUOwn - rhoUPhiOwn + pStarOwn*normal)
           /(rhoOwn*(SOwn - UvOwn)),
            Uf_
        );

        scalar f = (SOwn - UvOwn)/dS;
        phi = SStar*f;

        rhoUPhi =
            (SStar*(SOwn*rhoUOwn - rhoUPhiOwn) + SOwn*pStarOwn*normal)/dS;
        p = 0.5*(pStarNei + pStarOwn);
        rhoEPhi = SStar*(SOwn*rhoEOwn - rhoEPhiOwn + SOwn*pStarOwn)/dS;

        forAll(alphaPhis, phasei)
        {
            alphaPhis[phasei] = alphasOwn[phasei]*SStar*f;
            alphaRhoPhis[phasei] = alphasOwn[phasei]*rhosOwn[phasei]*SStar*f;
        }
    }
    else if (SNei > 0)
    {
        const scalar dS = SNei - SStar;

        this->save
        (
            facei,
            patchi,
            (SNei*rhoUNei - rhoUPhiNei + pStarNei*normal)
           /(rhoNei*(SNei - UvNei)),
            Uf_
        );

        scalar f = (SNei - UvNei)/dS;
        phi = SStar*f;

        rhoUPhi =
            (SStar*(SNei*rhoUNei - rhoUPhiNei) + SNei*pStarNei*normal)/dS;
        p = 0.5*(pStarNei + pStarOwn);
        rhoEPhi = SStar*(SNei*(rhoENei + pStarNei) - rhoEPhiNei)/dS;

        forAll(alphaPhis, phasei)
        {
            alphaPhis[phasei] = alphasNei[phasei]*SStar*f;
            alphaRhoPhis[phasei] = alphasNei[phasei]*rhosNei[phasei]*SStar*f;
        }
    }
    else
    {
        this->save(facei, patchi, UNei, Uf_);
        phi = UvNei;
        rhoUPhi = rhoUPhiNei;
        p = pNei;
        rhoEPhi = rhoEPhiNei;

        forAll(alphaPhis, phasei)
        {
            alphaPhis[phasei] = alphasNei[phasei]*UvNei;
            alphaRhoPhis[phasei] = alphasNei[phasei]*rhosNei[phasei]*UvNei;
        }
    }
    phi *= magSf;
    rhoUPhi *= magSf;
    rhoEPhi *= magSf;
    rhoEPhi += vMesh*magSf*p;
    forAll(alphaPhis, phasei)
    {
        alphaPhis[phasei] *= magSf;
        alphaRhoPhis[phasei] *= magSf;
    }
}


Foam::scalar Foam::fluxSchemes::HLLC::energyFlux
(
    const scalar& rhoOwn, const scalar& rhoNei,
    const vector& UOwn, const vector& UNei,
    const scalar& eOwn, const scalar& eNei,
    const scalar& pOwn, const scalar& pNei,
    const label facei, const label patchi
) const
{
    scalar SOwn = getValue(facei, patchi, SOwn_());
    scalar SNei = getValue(facei, patchi, SNei_());
    scalar SStar = getValue(facei, patchi, SStar_());
    scalar pStarOwn = getValue(facei, patchi, pStarOwn_());
    scalar pStarNei = getValue(facei, patchi, pStarNei_());
    scalar UvOwn = getValue(facei, patchi, UvOwn_());
    scalar UvNei = getValue(facei, patchi, UvNei_());
    scalar magSf = mag(getValue(facei, patchi, mesh_.Sf()));

    // Owner values
    const scalar rhoEOwn = rhoOwn*(eOwn + 0.5*magSqr(UOwn));
    const scalar rhoEPhiOwn = (rhoEOwn + pOwn)*UvOwn;

    // Neighbour values
    const scalar rhoENei = rhoNei*(eNei + 0.5*magSqr(UNei));
    const scalar rhoEPhiNei = (rhoENei + pNei)*UvNei;

    scalar rhoEPhi, p;
    if (SOwn > 0)
    {
        rhoEPhi = rhoEPhiOwn;
        p = pOwn;
    }
    else if (SStar > 0)
    {
        const scalar dS = SOwn - SStar;
        rhoEPhi = SStar*(SOwn*rhoEOwn - rhoEPhiOwn + SOwn*pStarOwn)/dS;
        p = 0.5*(pStarNei + pStarOwn);
    }
    else if (SNei > 0)
    {
        const scalar dS = SNei - SStar;
        rhoEPhi = SStar*(SNei*rhoENei - rhoEPhiNei + SNei*pStarNei)/dS;
        p = 0.5*(pStarNei + pStarOwn);
    }
    else
    {
        rhoEPhi = rhoEPhiNei;
        p = pNei;
    }

    return rhoEPhi*magSf + meshPhi(facei, patchi)*p;
}


Foam::scalar Foam::fluxSchemes::HLLC::interpolate
(
    const scalar& fOwn, const scalar& fNei,
    const label facei, const label patchi
) const
{
    scalar SOwn = getValue(facei, patchi, SOwn_());
    scalar SStar = getValue(facei, patchi, SStar_());

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
