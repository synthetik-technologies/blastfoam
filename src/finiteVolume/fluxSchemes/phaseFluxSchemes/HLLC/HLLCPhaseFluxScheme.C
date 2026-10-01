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

#include "HLLCPhaseFluxScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseFluxSchemes
{
    defineTypeNameAndDebug(HLLC, 0);
    addToRunTimeSelectionTable(phaseFluxScheme, HLLC, dictionary);
    addToRunTimeSelectionTable(phaseFluxScheme, HLLC, solid);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseFluxSchemes::HLLC::HLLC
(
    const surfaceScalarField& phi,
    const scalar residualAlpha
)
:
    phaseFluxScheme(phi, residualAlpha)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseFluxSchemes::HLLC::~HLLC()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::phaseFluxSchemes::HLLC::clear()
{
    phaseFluxScheme::clear();
    SOwn_.clear();
    SNei_.clear();
    SStar_.clear();
}

void Foam::phaseFluxSchemes::HLLC::createSavedFields()
{
    phaseFluxScheme::createSavedFields();
    if (SOwn_.valid())
    {
        SOwn_.ref() = Zero;
        SNei_.ref() = Zero;
        SStar_.ref() = Zero;
        return;
    }
    SOwn_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                fieldName("SOwn"),
                mesh_.time().name(),
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
                fieldName("SNei"),
                mesh_.time().name(),
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
                fieldName("SStar"),
                mesh_.time().name(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimVelocity, 0.0)
        )
    );
}


void Foam::phaseFluxSchemes::HLLC::calculateFluxes
(
    const scalar& alphaOwn, const scalar& alphaNei,
    const scalar& rhoO, const scalar& rhoN,
    const vector& UO, const vector& UN,
    const scalar& eO, const scalar& eN,
    const scalar& pO, const scalar& pN,
    const scalar& cO, const scalar& cN,
    const vector& Sf,
    scalar& phi,
    scalar& alphaRhoPhi,
    vector& alphaRhoUPhi,
    scalar& alphaRhoEPhi,
    const label facei, const label patchi
)
{
    scalar magSf = mag(Sf);
    vector normal = Sf/magSf;
    const bool ownValid = alphaOwn > 1e-6;
    const scalar rhoOwn = ownValid ? rhoO : 0.0;
    const vector UOwn = ownValid ? UO : vector::zero;
    const scalar eOwn = ownValid ? eO : 0.0;
    const scalar pOwn = ownValid ? pO : 0.0;
    const scalar cOwn = ownValid ? cO : 0.0;

    const bool neiValid = alphaNei > 1e-6;
    const scalar rhoNei = neiValid ? rhoN : 0.0;
    const vector UNei = neiValid ? UN : vector::zero;
    const scalar eNei = neiValid ? eN : 0.0;
    const scalar pNei = neiValid ? pN : 0.0;
    const scalar cNei = neiValid ? cN : 0.0;

    scalar EOwn = eOwn + 0.5*magSqr(UOwn);
    scalar ENei = eNei + 0.5*magSqr(UNei);

    const scalar phiMesh = meshPhi(facei, patchi);
    const scalar vMesh = phiMesh/magSf;
    scalar UvOwn((UOwn & normal) - vMesh);
    scalar UvNei((UNei & normal) - vMesh);

    // scalar SOwn(UvTilde - cTilde);
    // scalar SNei(UvTilde + cTilde);
    scalar SOwn(min(UvOwn - cOwn, UvNei - cNei));
    scalar SNei(max(UvOwn + cOwn, UvNei + cNei));

    const scalar denom = rhoOwn*(SOwn - UvOwn) - rhoNei*(SNei - UvNei);
    scalar SStar = 0.5*(UvOwn + UvNei);
    if (mag(denom) > 1e-8 && (cOwn > 1e-4 || cNei > 1e-4 || mag(pNei - pOwn) > 1e-4))
    {
        SStar =
            (
                pNei - pOwn
              + rhoOwn*UvOwn*(SOwn - UvOwn)
              - rhoNei*UvNei*(SNei - UvNei)
            )/denom;
    }
    scalar pStarOwn(pOwn + rhoOwn*(SOwn - UvOwn)*(SStar - UvOwn));
    scalar pStarNei(pNei + rhoNei*(SNei - UvNei)*(SStar - UvNei));

    this->save(facei, patchi, SOwn, SOwn_);
    this->save(facei, patchi, SNei, SNei_);
    this->save(facei, patchi, SStar, SStar_);

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

    scalar alpha;
    scalar rhoPhi;
    vector rhoUPhi;
    scalar rhoEPhi;
    scalar p;
    scalar f = 1.0;

    if (SOwn >= 0)
    {
        alpha = alphaOwn;
        phi = UvOwn;
        rhoPhi = rhoOwn*phi;
        rhoUPhi = rhoUPhiOwn;
        rhoEPhi = rhoEPhiOwn;
        p = pOwn;
        this->save(facei, patchi, UOwn, Uf_);
    }
    else if (SNei <= 0)
    {
        alpha = alphaNei;
        phi = UvNei;
        rhoPhi = rhoNei*phi;
        rhoUPhi = rhoUPhiNei;
        rhoEPhi = rhoEPhiNei;
        p = pNei;
        this->save(facei, patchi, UNei, Uf_);
    }
    else if (mag(SOwn - SStar) > small && SStar > 0)
    {
        f = (SOwn - UvOwn)/stabilise(SOwn - SStar, small);
        phi = SStar*f;
        alpha = alphaOwn;
        rhoPhi = rhoOwn*phi;

        vector rhoUStar = f*rhoOwn*(UOwn - (UvOwn - SStar)*normal);
        scalar rhoEStar =
            f*rhoOwn*EOwn + (pStarOwn*SStar - pOwn*UvOwn)/stabilise(SOwn - UvOwn, small);

        rhoUPhi = rhoUPhiOwn + SOwn*(rhoUStar - rhoUOwn);
        rhoEPhi = rhoEPhiOwn + SOwn*(rhoEStar - rhoEOwn);

        p = pStarOwn;
        this->save(facei, patchi, UOwn - (UvOwn - phi)*normal, Uf_);
    }
    else if (mag(SNei - SStar) > small && SStar < 0)
    {
        f = (SNei - UvNei)/stabilise(SNei - SStar, small);
        phi = SStar*f;
        alpha = alphaNei;
        rhoPhi = rhoNei*phi;

        vector rhoUStar = f*rhoNei*(UNei - (UvNei - SStar)*normal);
        scalar rhoEStar =
            f*rhoNei*ENei + (pStarNei*SStar - pNei*UvNei)/stabilise(SNei - UvNei, small);

        rhoUPhi = rhoUPhiNei + SNei*(rhoUStar - rhoUNei);
        rhoEPhi = rhoEPhiNei + SNei*(rhoEStar - rhoENei);

        p = pStarNei;
        this->save(facei, patchi, UNei - (UvNei - phi)*normal, Uf_);
    }
    else
    {
        alpha = 0.5*(alphaOwn + alphaNei);
        phi = 0.5*(UvOwn + UvNei);
        rhoPhi = 0.5*(rhoOwn*UvOwn + rhoNei*UvNei);
        rhoUPhi = 0.5*(rhoUPhiOwn + rhoUPhiNei);
        rhoEPhi = 0.5*(rhoEPhiOwn + rhoEPhiNei);
        p = 0.5*(pOwn + pNei);
        this->save(facei, patchi, 0.5*(UOwn + UNei), Uf_);
    }

    this->save(facei, patchi, alpha, alphaf_);
    this->save(facei, patchi, p, pf_);

    alphaRhoPhi = alpha*rhoPhi*magSf;
    alphaRhoUPhi = alpha*rhoUPhi*magSf;
    alphaRhoEPhi = alpha*rhoEPhi*magSf + phiMesh*alpha*p;
    phi *= magSf;
}


Foam::scalar Foam::phaseFluxSchemes::HLLC::calculateAlphaCorrector
(
    const scalar& alphaOwn, const scalar& alphaNei,
    const label facei, const label patchi
) const
{
    return 0.0;
    // NotImplemented;

    const scalar SOwn = this->getValue(facei, patchi, SOwn_);
    const scalar SNei = this->getValue(facei, patchi, SNei_);
    const scalar SStar = this->getValue(facei, patchi, SStar_);

    if (SOwn >= 0)
    {
        return 0.0;
    }
    else if (SStar > 0)
    {
        return
            SOwn
           *(
                (
                    SNei*alphaNei
                  - SOwn*alphaOwn
                  // - (SStar - SOwn)*(alphaNei - alphaOwn)
                )/(SNei - SOwn)
              - alphaOwn
            );
    }
    else if (SNei > 0)
    {
        return
            SNei
           *(
                (
                    SNei*alphaNei
                  - SOwn*alphaOwn
                  // - (SStar - SNei)*(alphaNei - alphaOwn)
                )/(SNei - SOwn)
              - alphaNei
            );
    }
    else
    {
        return 0.0;
    }
}


Foam::scalar Foam::phaseFluxSchemes::HLLC::interpolate
(
    const scalar& fOwn, const scalar& fNei,
    const label facei, const label patchi
) const
{
    return
            getValue(facei, patchi, SOwn_) > 0
         || getValue(facei, patchi, SStar_) > 0
          ? fOwn
          : fNei;
}


Foam::scalar Foam::phaseFluxSchemes::HLLC::calculateFlux
(
    const scalar& fOwn, const scalar& fNei,
    const scalar& phi,
    const label facei, const label patchi
) const
{
    return
        (
            getValue(facei, patchi, SOwn_) > 0
         || getValue(facei, patchi, SStar_) > 0
          ? fOwn
          : fNei
        )*phi;
}


// ************************************************************************* //
