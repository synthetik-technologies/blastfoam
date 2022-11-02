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

Foam::phaseFluxSchemes::HLLC::HLLC(const surfaceScalarField& phi)
:
    phaseFluxScheme(phi)
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
        return;
    }
    SOwn_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                fieldName("SOwn"),
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
                fieldName("SNei"),
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
                fieldName("SStar"),
                mesh_.time().timeName(),
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

    scalar SOwn(stabilise(min(UvOwn - cOwn, UvTilde - cTilde), small));
    scalar SNei(stabilise(max(UvNei + cNei, UvTilde + cTilde), small));

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

    // Owner values
    scalar alpha;
    scalar rho;
    vector U;
    scalar E;
    scalar p;
    scalar f = 1.0;

    if (SOwn >= 0)
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
        f = (SOwn - UvOwn)/(SOwn - SStar);
        alpha = alphaOwn;
        rho = rhoOwn*f;
        phi = SStar;
        U = (UOwn - UvOwn*normal) + SStar*normal;
        E = EOwn + (pStar*SStar - pOwn*UvOwn)/(rhoOwn*(SOwn - UvOwn));
        p = pStar;
    }
    else if (SNei > 0)
    {
        f = (SNei - UvNei)/(SNei - SStar);
        alpha = alphaNei;
        rho = rhoNei*f;
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
    alphaRhoEPhi = alphaPhi*(rho*E + p) + vMesh*magSf*p*alpha;
    phi *= f;
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
