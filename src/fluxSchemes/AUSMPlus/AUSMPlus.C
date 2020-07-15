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

#include "AUSMPlus.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluxSchemes
{
    defineTypeNameAndDebug(AUSMPlus, 0);
    addToRunTimeSelectionTable(fluxScheme, AUSMPlus, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxSchemes::AUSMPlus::AUSMPlus(const fvMesh& mesh)
:
    fluxScheme(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxSchemes::AUSMPlus::~AUSMPlus()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fluxSchemes::AUSMPlus::clear()
{
    fluxScheme::clear();
    phi_.clear();
}

void Foam::fluxSchemes::AUSMPlus::createSavedFields()
{
    fluxScheme::createSavedFields();
    if (phi_.valid())
    {
        return;
    }
    phi_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "AUSMPlus::phi",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimVelocity*dimArea, 0.0)
        )
    );
}


void Foam::fluxSchemes::AUSMPlus::calculateFluxes
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
    scalar HOwn(EOwn + pOwn/rhoOwn);

    scalar ENei = eNei + 0.5*magSqr(UNei);
    scalar HNei(ENei + pNei/rhoNei);

    const scalar vMesh(meshPhi(facei, patchi)/magSf);
    scalar UvOwn((UOwn & normal) - vMesh);
    scalar UvNei((UNei & normal) - vMesh);

    scalar c12(0.5*(cOwn + cNei));

    // Compute split Mach numbers
    scalar MaOwn(UvOwn/c12);
    scalar MaNei(UvNei/c12);
    scalar magMaOwn(mag(MaOwn));
    scalar magMaNei(mag(MaNei));

    scalar Ma4Own(max(MaOwn, 0.0));
    scalar P5Own(pos0(MaOwn));
    if (magMaOwn < 1)
    {
        Ma4Own = 0.25*sqr(MaOwn + 1.0) + beta_*sqr(sqr(MaOwn) - 1.0);
        P5Own =
            0.25*sqr(MaOwn + 1.0)*(2.0 - MaOwn)
          + alpha_*MaOwn*sqr(sqr(MaOwn) - 1.0);
    }

    scalar Ma4Nei(min(MaNei, 0.0));
    scalar P5Nei(neg(MaNei));
    if (magMaNei < 1)
    {
        Ma4Nei = -0.25*sqr(MaNei - 1.0) - beta_*sqr(sqr(MaNei) - 1.0);
        P5Nei =
            0.25*sqr(MaNei - 1.0)*(2.0 + MaNei)
          - alpha_*MaNei*sqr(sqr(MaNei) - 1.0);
    }
    scalar Ma12(Ma4Own + Ma4Nei);
    scalar P12(P5Own*pOwn + P5Nei*pNei);

    phi = magSf*c12*Ma12;

    this->save(facei, patchi, phi, phi_);

    scalar p;
    if (Ma12 >= 0)
    {
        this->save(facei, patchi, UOwn, Uf_);
        rhoPhi = rhoOwn;
        rhoUPhi = rhoOwn*UOwn;
        rhoEPhi = rhoOwn*HOwn;
        p = pOwn;
    }
    else
    {
        this->save(facei, patchi, UNei, Uf_);
        rhoPhi = rhoNei;
        rhoUPhi = rhoNei*UNei;
        rhoEPhi = rhoNei*HNei;
        p = pNei;
    }
    rhoPhi *= phi;
    rhoUPhi *= phi;
    rhoUPhi += P12*Sf;
    rhoEPhi *= phi;
    rhoEPhi += vMesh*magSf*p;
}


void Foam::fluxSchemes::AUSMPlus::calculateFluxes
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

    scalar c12(0.5*(cOwn + cNei));

    // Compute split Mach numbers
    scalar MaOwn(UvOwn/c12);
    scalar MaNei(UvNei/c12);
    scalar magMaOwn(mag(MaOwn));
    scalar magMaNei(mag(MaNei));

    scalar Ma4Own(max(MaOwn, 0.0));
    scalar P5Own(pos0(MaOwn));
    if (magMaOwn < 1)
    {
        Ma4Own = 0.25*sqr(MaOwn + 1.0) + beta_*sqr(sqr(MaOwn) - 1.0);
        P5Own =
            0.25*sqr(MaOwn + 1.0)*(2.0 - MaOwn)
          + alpha_*MaOwn*sqr(sqr(MaOwn) - 1.0);
    }

    scalar Ma4Nei(min(MaNei, 0.0));
    scalar P5Nei(neg(MaNei));
    if (magMaNei < 1)
    {
        Ma4Nei = -0.25*sqr(MaNei - 1.0) - beta_*sqr(sqr(MaNei) - 1.0);
        P5Nei =
            0.25*sqr(MaNei - 1.0)*(2.0 + MaNei)
          - alpha_*MaNei*sqr(sqr(MaNei) - 1.0);
    }
    scalar Ma12(Ma4Own + Ma4Nei);
    scalar P12(P5Own*pOwn + P5Nei*pNei);

    phi = magSf*c12*Ma12;

    this->save(facei, patchi, phi, phi_);

    scalar p;
    if (Ma12 >= 0)
    {
        this->save(facei, patchi, UOwn, Uf_);
        rhoUPhi = rhoOwn*UOwn;
        rhoEPhi = rhoOwn*EOwn + pOwn;
        p = pOwn;
        forAll(alphasOwn, phasei)
        {
            alphaPhis[phasei] = alphasOwn[phasei];
            alphaRhoPhis[phasei] = alphasOwn[phasei]*rhosOwn[phasei];
        }
    }
    else
    {
        this->save(facei, patchi, UNei, Uf_);
        rhoUPhi = rhoNei*UNei;
        rhoEPhi = rhoNei*ENei + pNei;
        p = pNei;
        forAll(alphasNei, phasei)
        {
            alphaPhis[phasei] = alphasNei[phasei];
            alphaRhoPhis[phasei] = alphasNei[phasei]*rhosNei[phasei];
        }
    }
    rhoUPhi *= phi;
    rhoUPhi += P12*Sf;
    rhoEPhi *= phi;
    rhoEPhi += vMesh*magSf*p;

    forAll(alphasOwn, phasei)
    {
        alphaPhis[phasei] *= phi;
        alphaRhoPhis[phasei] *= phi;
    }
}


Foam::scalar Foam::fluxSchemes::AUSMPlus::energyFlux
(
    const scalar& rhoOwn, const scalar& rhoNei,
    const vector& UOwn, const vector& UNei,
    const scalar& eOwn, const scalar& eNei,
    const scalar& pOwn, const scalar& pNei,
    const label facei, const label patchi
) const
{
    scalar phi = getValue(facei, patchi, phi_());
    if ( phi >= 0)
    {
        return
            phi*(rhoOwn*(eOwn + 0.5*magSqr(UOwn)) + pOwn)
          + meshPhi(facei, patchi)*pOwn;
    }
    else
    {
        return
            phi*(rhoNei*(eNei + 0.5*magSqr(UNei)) + pNei)
          + meshPhi(facei, patchi)*pNei;
    }
}


Foam::scalar Foam::fluxSchemes::AUSMPlus::interpolate
(
    const scalar& fOwn, const scalar& fNei,
    const label facei, const label patchi
) const
{

    if (getValue(facei, patchi, phi_()) >= 0)
    {
        return fOwn;
    }
    else
    {
        return fNei;
    }
}
