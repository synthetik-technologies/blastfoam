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

#include "HLL.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluxSchemes
{
    defineTypeNameAndDebug(HLL, 0);
    addToRunTimeSelectionTable(fluxScheme, HLL, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxSchemes::HLL::HLL
(
    const fvMesh& mesh
)
:
    fluxScheme(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxSchemes::HLL::~HLL()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluxSchemes::HLL::clear()
{
    fluxScheme::clear();
    if (SOwn_.valid())
    {
        return;
    }

    SOwn_.clear();
    SNei_.clear();
    UvOwn_.clear();
    UvNei_.clear();
}

void Foam::fluxSchemes::HLL::createSavedFields()
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
                "HLL::SOwn",
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
                "HLL::SNei",
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
                "HLL::UvOwn",
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
                "HLL::UvNei",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimVelocity, 0.0)
        )
    );
}


void Foam::fluxSchemes::HLL::calculateFluxes
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

    scalar wOwn(sqrt(rhoOwn)/(sqrt(rhoOwn) + sqrt(rhoNei)));
    scalar wNei(1.0 - wOwn);

    scalar cTilde(cOwn*wOwn + cNei*wNei);
    scalar UvTilde(UvOwn*wOwn + UvNei*wNei);

    scalar SOwn(min(UvOwn - cOwn, UvTilde - cTilde));
    scalar SNei(max(UvNei + cNei, UvTilde + cTilde));

    this->save(facei, patchi, SOwn, SOwn_);
    this->save(facei, patchi, SNei, SNei_);
    this->save(facei, patchi, UvOwn, UvOwn_);
    this->save(facei, patchi, UvNei, UvNei_);

    scalar p;
    if (SOwn >= 0)
    {
        phi = this->save(facei, patchi, UOwn, Uf_) & normal;
        rhoPhi = rhoOwn*UvOwn;
        rhoUPhi = UOwn*rhoOwn*UvOwn + pOwn*normal;
        rhoEPhi = UvOwn*rhoOwn*HOwn;
        p = pOwn;
    }
    else if (SOwn < 0 && SNei >= 0)
    {
        vector rhoUOwn = rhoOwn*UOwn;
        vector rhoUNei = rhoNei*UNei;
        scalar rhoPhiOwn = rhoOwn*UvOwn;
        scalar rhoPhiNei = rhoNei*UvNei;
        vector rhoUPhiOwn = rhoUOwn*UvOwn + pOwn*normal;
        vector rhoUPhiNei = rhoUNei*UvNei + pNei*normal;
        phi =
            this->save
            (
                facei,
                patchi,
                (SNei*rhoUNei - SOwn*rhoUOwn + rhoUPhiOwn - rhoUPhiNei)
               /(SNei*rhoNei - SOwn*rhoOwn + rhoPhiOwn - rhoPhiNei),
                Uf_
            ) & normal;

        rhoPhi =
            (
                SNei*rhoOwn*UvOwn - SOwn*rhoNei*UvNei
              + SOwn*SNei*(rhoNei - rhoOwn)
            )/(SNei - SOwn);

        rhoUPhi =
            (
                SNei*rhoUPhiOwn - SOwn*rhoUPhiNei
              + SOwn*SNei*(rhoUNei - rhoUOwn)
            )/(SNei - SOwn);

        rhoEPhi =
            (
                SNei*rhoOwn*HOwn*UvOwn
              - SOwn*rhoNei*HNei*UvNei
              + SOwn*SNei
               *(
                    rhoNei*ENei
                  - rhoOwn*EOwn
                )
            )/(SNei - SOwn);
        p = (SNei*pOwn - SOwn*pNei)/(SNei - SOwn);
    }
    else
    {
        phi = this->save(facei, patchi, UNei, Uf_) & normal;
        rhoPhi = rhoNei*UvNei;
        rhoUPhi = UNei*rhoNei*UvNei + pNei*normal;
        rhoEPhi = rhoNei*HNei*UvNei;
        p = pNei;
    }
    phi *= magSf;
    rhoPhi *= magSf;
    rhoUPhi *= magSf;
    rhoEPhi *= magSf;
    rhoEPhi += vMesh*magSf*p;
}


void Foam::fluxSchemes::HLL::calculateFluxes
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
    scalar HOwn(EOwn + pOwn/rhoOwn);

    scalar ENei = eNei + 0.5*magSqr(UNei);
    scalar HNei(ENei + pNei/rhoNei);

    const scalar vMesh(meshPhi(facei, patchi)/magSf);
    scalar UvOwn((UOwn & normal) - vMesh);
    scalar UvNei((UNei & normal) - vMesh);

    scalar wOwn(sqrt(rhoOwn)/(sqrt(rhoOwn) + sqrt(rhoNei)));
    scalar wNei(1.0 - wOwn);

    scalar cTilde(cOwn*wOwn + cNei*wNei);
    scalar UvTilde(UvOwn*wOwn + UvNei*wNei);

    scalar SOwn(min(UvOwn - cOwn, UvTilde - cTilde));
    scalar SNei(max(UvNei + cNei, UvTilde + cTilde));

    this->save(facei, patchi, SOwn, SOwn_);
    this->save(facei, patchi, SNei, SNei_);
    this->save(facei, patchi, UvOwn, UvOwn_);
    this->save(facei, patchi, UvNei, UvNei_);

    scalar p;
    if (SOwn >= 0)
    {
        phi = this->save(facei, patchi, UOwn, Uf_) & normal;
        rhoUPhi = UOwn*rhoOwn*UvOwn + pOwn*normal;
        rhoEPhi = UvOwn*rhoOwn*HOwn;
        p = pOwn;

        forAll(alphasOwn, phasei)
        {
            alphaPhis[phasei] = alphasOwn[phasei]*UvOwn;
            alphaRhoPhis[phasei] = alphasOwn[phasei]*rhosOwn[phasei]*UvOwn;
        }
    }
    else if (SOwn < 0 && SNei >= 0)
    {
        vector rhoUOwn = rhoOwn*UOwn;
        vector rhoUNei = rhoNei*UNei;
        scalar rhoPhiOwn = rhoOwn*UvOwn;
        scalar rhoPhiNei = rhoNei*UvNei;
        vector rhoUPhiOwn = rhoUOwn*UvOwn + pOwn*normal;
        vector rhoUPhiNei = rhoUNei*UvNei + pNei*normal;
        phi =
            this->save
            (
                facei,
                patchi,
                (SNei*rhoUNei - SOwn*rhoUOwn + rhoUPhiOwn - rhoUPhiNei)
               /(SNei*rhoNei - SOwn*rhoOwn + rhoPhiOwn - rhoPhiNei),
                Uf_
            ) & normal;

        rhoUPhi =
            (
                SNei*(rhoUOwn*UvOwn + pOwn*normal)
              - SOwn*(rhoUNei*UvNei + pNei*normal)
              + SOwn*SNei*(rhoUNei - rhoUOwn)
            )/(SNei - SOwn);
        p = (SNei*pOwn - SOwn*pNei)/(SNei - SOwn);

        rhoEPhi =
            (
                SNei*rhoOwn*HOwn*UvOwn
              - SOwn*rhoNei*HNei*UvNei
              + SOwn*SNei
               *(
                    rhoNei*ENei
                  - rhoOwn*EOwn
                )
            )/(SNei - SOwn);

        forAll(alphasOwn, phasei)
        {
            alphaPhis[phasei] =
            (
                SNei*alphasOwn[phasei]*UvOwn - SOwn*alphasNei[phasei]*UvNei
              + SOwn*SNei*(alphasNei[phasei] - alphasOwn[phasei])
            )/(SNei - SOwn);

            scalar alphaRhoOwn = alphasOwn[phasei]*rhosOwn[phasei];
            scalar alphaRhoNei = alphasNei[phasei]*rhosNei[phasei];
            alphaRhoPhis[phasei] =
                (
                    SNei*alphaRhoOwn*UvOwn - SOwn*alphaRhoNei*UvNei
                  + SOwn*SNei*(alphaRhoNei - alphaRhoOwn)
                )/(SNei - SOwn);
        }
    }
    else
    {
        phi = this->save(facei, patchi, UNei, Uf_) & normal;
        rhoUPhi = UNei*rhoNei*UvNei + pNei*normal;
        rhoEPhi = rhoNei*HNei*UvNei;
        p = pNei;

        forAll(alphasOwn, phasei)
        {
            alphaPhis[phasei] = alphasNei[phasei]*UvNei;
            alphaRhoPhis[phasei] = alphasNei[phasei]*rhosNei[phasei]*UvNei;
        }
    }
    phi *= magSf;
    rhoUPhi *= magSf;
    rhoEPhi *= magSf;
    rhoEPhi += vMesh*magSf*p;
    forAll(alphasOwn, phasei)
    {
        alphaPhis[phasei] *= magSf;
        alphaRhoPhis[phasei] *= magSf;
    }
}


Foam::scalar Foam::fluxSchemes::HLL::energyFlux
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
    scalar UvOwn = getValue(facei, patchi, UvOwn_());
    scalar UvNei = getValue(facei, patchi, UvNei_());
    scalar magSf = mag(getValue(facei, patchi, mesh_.Sf()));

    scalar EOwn = eOwn + 0.5*magSqr(UOwn);
    scalar HOwn = EOwn + pOwn/rhoOwn;

    scalar ENei = eNei + 0.5*magSqr(UNei);
    scalar HNei = ENei + pNei/rhoNei;

    scalar phi, p;
    if (SOwn >= 0)
    {
        phi = UvOwn*rhoOwn*HOwn + meshPhi(facei, patchi)*pOwn;
        p = pOwn;
    }
    else if (SOwn < 0 && SNei >= 0)
    {
        phi =
            (
                SNei*rhoOwn*HOwn*UvOwn
              - SOwn*rhoNei*HNei*UvNei
              + SOwn*SNei
               *(
                    rhoNei*ENei
                  - rhoOwn*EOwn
                )
            )/(SNei - SOwn);
        p = (SNei*pOwn - SOwn*pNei)/(SNei - SOwn);
    }
    else
    {
        phi = rhoNei*HNei*UvNei + meshPhi(facei, patchi)*pNei;
        p = pNei;
    }
    return phi*magSf + meshPhi(facei, patchi)*p;
}


Foam::scalar Foam::fluxSchemes::HLL::interpolate
(
    const scalar& fOwn, const scalar& fNei,
    const label facei, const label patchi
) const
{
    scalar SOwn = getValue(facei, patchi, SOwn_());
    scalar SNei = getValue(facei, patchi, SNei_());
    scalar UvOwn = getValue(facei, patchi, UvOwn_());
    scalar UvNei = getValue(facei, patchi, UvNei_());

    if (SOwn >= 0)
    {
        return fOwn;
    }
    else if (SOwn < 0 && SNei >= 0)
    {
        return
            (SNei*fNei - SOwn*fOwn + fOwn*UvOwn + fNei*UvNei)/(SNei - SOwn);
    }
    else
    {
        return fNei;
    }
}

// ************************************************************************* //
