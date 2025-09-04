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

#include "HLLPhaseFluxScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseFluxSchemes
{
    defineTypeNameAndDebug(HLL, 0);
    addToRunTimeSelectionTable(phaseFluxScheme, HLL, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseFluxSchemes::HLL::HLL
(
    const surfaceScalarField& phi,
    const scalar residualAlpha
)
:
    phaseFluxScheme(phi, residualAlpha)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseFluxSchemes::HLL::~HLL()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::phaseFluxSchemes::HLL::clear()
{
    phaseFluxScheme::clear();
    SOwn_.clear();
    SNei_.clear();
    UvOwn_.clear();
    UvNei_.clear();
}

void Foam::phaseFluxSchemes::HLL::createSavedFields()
{
    phaseFluxScheme::createSavedFields();
    if (SOwn_.valid())
    {
        SOwn_.ref() = Zero;
        SNei_.ref() = Zero;
        UvOwn_.ref() = Zero;
        UvNei_.ref() = Zero;
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
    UvOwn_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                fieldName("UvOwn"),
                mesh_.time().name(),
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
                fieldName("UvNei"),
                mesh_.time().name(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimVelocity, 0.0)
        )
    );
}


void Foam::phaseFluxSchemes::HLL::calculateFluxes
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
    const scalar rhoOwn = ownValid ? rhoO : small;
    const vector UOwn = ownValid ? UO : vector::zero;
    const scalar eOwn = ownValid ? eO : 0.0;
    const scalar pOwn = ownValid ? pO : 0.0;
    const scalar cOwn = ownValid ? cO : 0.0;

    const bool neiValid = alphaNei > 1e-6;
    const scalar rhoNei = neiValid ? rhoN : small;
    const vector UNei = neiValid ? UN : vector::zero;
    const scalar eNei = neiValid ? eN : 0.0;
    const scalar pNei = neiValid ? pN : 0.0;
    const scalar cNei = neiValid ? cN : 0.0;

    const scalar phiMesh = meshPhi(facei, patchi);
    const scalar vMesh = phiMesh/magSf;
    scalar UvOwn((UOwn & normal) - vMesh);
    scalar UvNei((UNei & normal) - vMesh);

    scalar EOwn = eOwn + 0.5*magSqr(UOwn);
    scalar HOwn(EOwn + pOwn/rhoOwn);

    scalar ENei = eNei + 0.5*magSqr(UNei);
    scalar HNei(ENei + pNei/rhoNei);

    scalar SOwn(min(UvOwn - cOwn, UvNei - cNei));
    scalar SNei(max(UvOwn + cOwn, UvNei + cNei));

    this->save(facei, patchi, SOwn, SOwn_);
    this->save(facei, patchi, SNei, SNei_);
    this->save(facei, patchi, UvOwn, UvOwn_);
    this->save(facei, patchi, UvNei, UvNei_);

    // Owner values
    scalar alpha, p;
    vector U;

    if (SOwn >= 0)
    {
        alphaRhoPhi = alphaOwn*UvOwn*rhoOwn*magSf;
        alphaRhoUPhi = alphaRhoPhi*UOwn + alphaOwn*pOwn*Sf;
        alphaRhoEPhi = alphaRhoPhi*HOwn + phiMesh*alphaOwn*pOwn;

        alpha = alphaOwn;
        U = UOwn;
        p = pOwn;
    }
    else if (SOwn < 0 && SNei > 0)
    {
        const scalar rDeltaS(1.0/(SNei - SOwn));

        const scalar alphaRhoOwn = alphaOwn*rhoOwn;
        const scalar alphaRhoNei = alphaNei*rhoNei;

        const vector alphaRhoUOwn = alphaRhoOwn*UOwn;
        const vector alphaRhoUNei = alphaRhoNei*UNei;

        const scalar alphaRhoEOwn = alphaRhoOwn*EOwn;
        const scalar alphaRhoENei = alphaRhoNei*ENei;

        const scalar alphaRhoPhiOwn = alphaRhoOwn*UvOwn*magSf;
        const scalar alphaRhoPhiNei = alphaRhoNei*UvNei*magSf;

        const vector alphaRhoUPhiOwn = alphaRhoPhiOwn*UOwn + alphaOwn*pOwn*Sf;
        const vector alphaRhoUPhiNei = alphaRhoPhiNei*UNei + alphaNei*pNei*Sf;

        const scalar alphaRhoEPhiOwn =
            alphaRhoPhiOwn*HOwn + phiMesh*alphaOwn*pOwn;
        const scalar alphaRhoEPhiNei =
            alphaRhoPhiNei*HNei + phiMesh*alphaNei*pNei;

        alphaRhoPhi =
            (
                SNei*alphaRhoPhiOwn - SOwn*alphaRhoPhiNei
              + SOwn*SNei*(alphaRhoNei - alphaRhoOwn)*magSf
            )*rDeltaS;

        alphaRhoUPhi =
            (
                SNei*alphaRhoUPhiOwn - SOwn*alphaRhoUPhiNei
              + SOwn*SNei*(alphaRhoUNei - alphaRhoUOwn)*magSf
            )*rDeltaS;

        alphaRhoEPhi =
            (
                SNei*alphaRhoEPhiOwn - SOwn*alphaRhoEPhiNei
              + SOwn*SNei*(alphaRhoENei - alphaRhoEOwn)*magSf
            )*rDeltaS;

        alpha = (SNei*alphaOwn - SOwn*alphaNei)*rDeltaS;

        U = // alphaRhoU_hll / alphaRho_hll
            (
                (SNei*alphaRhoUNei - SOwn*alphaRhoUOwn)*magSf
              + alphaRhoUPhiOwn - alphaRhoUPhiNei
            )
            /(
                (SNei*alphaRhoNei - SOwn*alphaRhoOwn)*magSf
              + alphaRhoPhiOwn - alphaRhoPhiNei
            );
        p = 0.5*(pOwn + pNei);
    }
    else
    {
        alphaRhoPhi = alphaNei*UvNei*rhoNei*magSf;
        alphaRhoUPhi = alphaRhoPhi*UNei + alphaNei*pNei*Sf;
        alphaRhoEPhi = alphaRhoPhi*HNei + phiMesh*alphaNei*pNei;

        alpha = alphaNei;
        U = UNei;
        p = pNei;
    }

    this->save(facei, patchi, alpha, alphaf_);
    phi = this->save(facei, patchi, U, Uf_) & Sf;
    this->save(facei, patchi, p, pf_);
}

Foam::scalar Foam::phaseFluxSchemes::HLL::calculateAlphaCorrector
(
    const scalar& alphaOwn, const scalar& alphaNei,
    const label facei, const label patchi
) const
{
    const scalar SOwn = this->getValue(facei, patchi, SOwn_);
    const scalar SNei = this->getValue(facei, patchi, SNei_);

    if (SOwn >= 0)
    {
        return 0.0;
    }
    else if (SOwn < 0 && SNei > 0)
    {
        return SOwn*SNei*(alphaNei - alphaOwn)/(SNei - SOwn);
    }
    else
    {
        return 0.0;
    }
}


Foam::scalar Foam::phaseFluxSchemes::HLL::interpolate
(
    const scalar& fOwn, const scalar& fNei,
    const label facei, const label patchi
) const
{
    scalar SOwn = getValue(facei, patchi, SOwn_());
    scalar SNei = getValue(facei, patchi, SNei_());

    if (SOwn >= 0)
    {
        return fOwn;
    }
    else if (SOwn < 0 && SNei >= 0)
    {
        scalar UvOwn = getValue(facei, patchi, UvOwn_());
        scalar UvNei = getValue(facei, patchi, UvNei_());
        return
            (SNei*fNei - SOwn*fOwn + fOwn*UvOwn - fNei*UvNei)/(SNei - SOwn);
    }
    else
    {
        return fNei;
    }
}


Foam::scalar Foam::phaseFluxSchemes::HLL::calculateFlux
(
    const scalar& fOwn, const scalar& fNei,
    const scalar& phi,
    const label facei, const label patchi
) const
{
    scalar SOwn = getValue(facei, patchi, SOwn_);
    scalar SNei = getValue(facei, patchi, SNei_);
    if (SOwn >= 0)
    {
        return fOwn*phi;
    }
    else if (SOwn < 0 && SNei >= 0)
    {
        return
            (
                phi*(SNei*fOwn - SOwn*fNei)
              + SOwn*SNei*(fNei - fOwn)*getValue(facei, patchi, mesh_.magSf())
            )/(SNei - SOwn);
    }
    else
    {
        return fNei*phi;
    }
}


// ************************************************************************* //
