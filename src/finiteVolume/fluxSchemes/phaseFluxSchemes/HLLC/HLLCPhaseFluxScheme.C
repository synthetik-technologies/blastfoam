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


    scalar SStar
    (
        (
            pNei - pOwn
          + rhoOwn*UvOwn*(SOwn - UvOwn)
          - rhoNei*UvNei*(SNei - UvNei)
        )
       /stabilise(rhoOwn*(SOwn - UvOwn) - rhoNei*(SNei - UvNei), small)
    );
    scalar pStar =
        0.5
       *(
            pOwn + rhoOwn*(SOwn - UvOwn)*(SStar - UvOwn)
          + pNei + rhoNei*(SNei - UvNei)*(SStar - UvNei)
        );

    this->save(facei, patchi, SOwn, SOwn_);
    this->save(facei, patchi, SNei, SNei_);
    this->save(facei, patchi, SStar, SStar_);

    // Owner values
    scalar alpha, rho, E, p;
    vector U;
    scalar f = 1.0;

    if (SOwn >= 0)
    {
        alpha = alphaOwn;
        rho = rhoOwn;
        phi = UvOwn*magSf;
        U = UOwn;
        E = EOwn;
        p = pOwn;
    }
    else if (SStar > 0)
    {
        const scalar rDeltaS = 1.0/(SOwn - SStar);
        f = (SOwn - UvOwn)*rDeltaS;

        alpha = alphaOwn;
        rho = rhoOwn*f;
        phi = SStar*magSf;
        U = (UOwn - UvOwn*normal) + SStar*normal;
        E =
            EOwn
          + (SStar - UvOwn)
           *(SStar + pOwn/stabilise(rhoOwn*(SOwn - UvOwn), small));
          // + (pStar*SStar - pOwn*UvOwn)
           // /stabilise(rhoOwn*(SOwn - UvOwn), small);
        p = pStar;
    }
    else if (SNei > 0)
    {
        const scalar rDeltaS = 1.0/(SNei - SStar);
        f = (SNei - UvNei)*rDeltaS;
        alpha = alphaNei;
        rho = rhoNei*f;
        phi = SStar*magSf;
        U = (UNei - UvNei*normal) + SStar*normal;
        E =
            ENei
          + (SStar - UvNei)
           *(SStar + pNei/stabilise(rhoNei*(SNei - UvNei), small));
          // + (pStar*SStar - pNei*UvNei)
           // /stabilise(rhoNei*(SNei - UvNei), small);
        p = pStar;
    }
    else
    {
        alpha = alphaNei;
        rho = rhoNei;
        phi = UvNei*magSf;
        U = UNei;
        E = ENei;
        p = pNei;
    }

    this->save(facei, patchi, alpha, alphaf_);
    this->save(facei, patchi, U, Uf_);
    this->save(facei, patchi, p, pf_);

    alphaRhoPhi = alpha*rho*phi;
    alphaRhoUPhi = alphaRhoPhi*U + alpha*p*Sf;
    alphaRhoEPhi = alpha*phi*(rho*E + p) + phiMesh*alpha*p;
    phi *= f;
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
