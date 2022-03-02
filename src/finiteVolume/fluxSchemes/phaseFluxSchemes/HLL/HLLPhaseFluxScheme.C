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
    const fvMesh& mesh,
    const word& name
)
:
    phaseFluxScheme(mesh, name)
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
        return;
    }
    SOwn_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                IOobject::groupName("HLL::SOwn", this->group()),
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
                IOobject::groupName("HLL::SNei", this->group()),
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
                IOobject::groupName("HLL::UvOwn", this->group()),
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
                IOobject::groupName("HLL::UvNei", this->group()),
                mesh_.time().timeName(),
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

    const scalar vMesh(meshPhi(facei, patchi)/magSf);
    scalar UvOwn((UOwn & normal) - vMesh);
    scalar UvNei((UNei & normal) - vMesh);

    scalar rhoOwn = max(rhoO, 1e-10);
    scalar rhoNei = max(rhoN, 1e-10);

    scalar EOwn = eOwn + 0.5*magSqr(UOwn);
    scalar HOwn(EOwn + pOwn/rhoOwn);

    scalar ENei = eNei + 0.5*magSqr(UNei);
    scalar HNei(ENei + pNei/rhoNei);

    scalar SOwn(stabilise(min(UvOwn - cOwn, UvNei - cNei), small));
    scalar SNei(stabilise(max(UvOwn + cOwn, UvNei + cNei), small));

    this->save(facei, patchi, SOwn, SOwn_);
    this->save(facei, patchi, SNei, SNei_);
    this->save(facei, patchi, UvOwn, UvOwn_);
    this->save(facei, patchi, UvNei, UvNei_);

    // Owner values
    scalar alpha;
    vector U;
    scalar p;

    if (SOwn >= 0)
    {
        alphaPhi = alphaOwn*UvOwn;
        alphaRhoPhi = alphaPhi*rhoOwn;
        alphaRhoUPhi = alphaRhoPhi*UOwn + alphaOwn*pOwn*normal;
        alphaRhoEPhi = alphaRhoPhi*HOwn;

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

        const scalar alphaPhiOwn = alphaOwn*UvOwn;
        const scalar alphaPhiNei = alphaNei*UvNei;

        const scalar alphaRhoPhiOwn = alphaRhoOwn*UvOwn;
        const scalar alphaRhoPhiNei = alphaRhoNei*UvNei;

        const vector alphaRhoUPhiOwn = alphaRhoUOwn*UvOwn + alphaOwn*pOwn*normal;
        const vector alphaRhoUPhiNei = alphaRhoUNei*UvNei + alphaNei*pNei*normal;

        const scalar alphaRhoEPhiOwn = alphaRhoPhiOwn*HOwn;
        const scalar alphaRhoEPhiNei = alphaRhoPhiNei*HNei;

        alphaPhi =
            (
                SNei*alphaPhiOwn - SOwn*alphaPhiNei
              + SOwn*SNei*(alphaNei - alphaOwn)
            )*rDeltaS;

        alphaRhoPhi =
            (
                SNei*alphaRhoPhiOwn - SOwn*alphaRhoPhiNei
              + SOwn*SNei*(alphaRhoNei - alphaRhoOwn)
            )*rDeltaS;

        alphaRhoUPhi =
            (
                SNei*alphaRhoUPhiOwn - SOwn*alphaRhoUPhiNei
              + SOwn*SNei*(alphaRhoUNei - alphaRhoUOwn)
            )*rDeltaS;

        alphaRhoEPhi =
            (
                SNei*alphaRhoEPhiOwn - SOwn*alphaRhoEPhiNei
              + SOwn*SNei*(alphaRhoENei - alphaRhoEOwn)
            )*rDeltaS;

        alpha = 0.5*(alphaOwn + alphaNei);
        U = 0.5*(UOwn + UNei);
//             (
//                 SNei*alphaRhoUNei - SOwn*alphaRhoUOwn
//               + alphaRhoUPhiOwn - alphaRhoUPhiNei
//             )
//             /(
//                 SNei*alphaRhoNei - SOwn*alphaRhoOwn
//               + alphaRhoPhiOwn - alphaRhoPhiNei
//             );
        p = 0.5*(pOwn + pNei);
    }
    else
    {
        phi = this->save(facei, patchi, UNei, Uf_) & normal;
        alphaPhi = alphaNei*UvNei;
        alphaRhoPhi = alphaPhi*rhoNei;
        alphaRhoUPhi = alphaRhoPhi*UNei + alphaNei*pNei*normal;
        alphaRhoEPhi = alphaRhoPhi*HNei;

        alpha = alphaNei;
        U = UNei;
        p = pNei;
    }

    this->save(facei, patchi, alpha, alphaf_);
    phi = this->save(facei, patchi, U, Uf_) & normal;
    this->save(facei, patchi, p, pf_);

    phi *= magSf;
    alphaPhi *= magSf;
    alphaRhoPhi *= magSf;
    alphaRhoUPhi *= magSf;
    alphaRhoEPhi *= magSf;
    alphaRhoEPhi += vMesh*magSf*alpha*p;
}


void Foam::phaseFluxSchemes::HLL::calculateFluxes
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

    const scalar vMesh(meshPhi(facei, patchi)/magSf);
    scalar UvOwn((UOwn & normal) - vMesh);
    scalar UvNei((UNei & normal) - vMesh);

    scalar rhoOwn = max(rhoO, 1e-10);
    scalar rhoNei = max(rhoN, 1e-10);

    scalar EOwn = eOwn + 0.5*magSqr(UOwn);
    scalar HOwn(EOwn + pOwn/rhoOwn);

    scalar ENei = eNei + 0.5*magSqr(UNei);
    scalar HNei(ENei + pNei/rhoNei);

    scalar SOwn(stabilise(min(UvOwn - cOwn, UvNei - cNei), small));
    scalar SNei(stabilise(max(UvOwn + cOwn, UvNei + cNei), small));

    this->save(facei, patchi, SOwn, SOwn_);
    this->save(facei, patchi, SNei, SNei_);
    this->save(facei, patchi, UvOwn, UvOwn_);
    this->save(facei, patchi, UvNei, UvNei_);

    // Owner values
    scalar alpha;
    vector U;
    scalar p;

    if (SOwn >= 0)
    {
        phi = this->save(facei, patchi, UOwn, Uf_) & normal;
        scalar alphaRhoPhi = alphaOwn*rhoOwn*UvOwn;
        alphaRhoUPhi = alphaRhoPhi*UOwn + alphaOwn*pOwn*normal;
        alphaRhoEPhi = alphaRhoPhi*HOwn;

        forAll(alphasOwn, phasei)
        {
            alphaPhis[phasei] = alphasOwn[phasei]*UvOwn;
            alphaRhoPhis[phasei] = alphaPhis[phasei]*rhosOwn[phasei];
        }

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

        const scalar alphaRhoPhiOwn = alphaRhoOwn*UvOwn;
        const scalar alphaRhoPhiNei = alphaRhoNei*UvNei;

        const vector alphaRhoUPhiOwn = alphaRhoUOwn*UvOwn + alphaOwn*pOwn*normal;
        const vector alphaRhoUPhiNei = alphaRhoUNei*UvNei + alphaNei*pNei*normal;

        const scalar alphaRhoEPhiOwn = alphaRhoPhiOwn*HOwn;
        const scalar alphaRhoEPhiNei = alphaRhoPhiNei*HNei;

        alphaRhoUPhi =
            (
                SNei*alphaRhoUPhiOwn - SOwn*alphaRhoUPhiNei
              + SOwn*SNei*(alphaRhoUNei - alphaRhoUOwn)
            )*rDeltaS;

        alphaRhoEPhi =
            (
                SNei*alphaRhoEPhiOwn - SOwn*alphaRhoEPhiNei
              + SOwn*SNei*(alphaRhoENei - alphaRhoEOwn)
            )*rDeltaS;

        forAll(alphasOwn, phasei)
        {
            const scalar alphaRhoIOwn = alphasOwn[phasei]*rhosOwn[phasei];
            const scalar alphaRhoINei = alphasNei[phasei]*rhosNei[phasei];

            const scalar alphaPhiIOwn = alphasOwn[phasei]*UvOwn;
            const scalar alphaPhiINei = alphasNei[phasei]*UvNei;

            const scalar alphaRhoPhiIOwn = alphaRhoOwn*UvOwn;
            const scalar alphaRhoPhiINei = alphaRhoNei*UvNei;

            alphaPhis[phasei] =
                (
                    SNei*alphaPhiIOwn - SOwn*alphaPhiINei
                  + SOwn*SNei*(alphasNei[phasei] - alphasOwn[phasei])
                )*rDeltaS;

            alphaRhoPhis[phasei] =
                (
                    SNei*alphaRhoPhiIOwn - SOwn*alphaRhoPhiINei
                  + SOwn*SNei*(alphaRhoINei - alphaRhoIOwn)
                )*rDeltaS;
        }

        alpha = 0.5*(alphaOwn + alphaNei);
        U = 0.5*(UOwn + UNei);
//             (
//                 SNei*alphaRhoUNei - SOwn*alphaRhoUOwn
//               + alphaRhoUPhiOwn - alphaRhoUPhiNei
//             )
//             /(
//                 SNei*alphaRhoNei - SOwn*alphaRhoOwn
//               + alphaRhoPhiOwn - alphaRhoPhiNei
//             );
        p = 0.5*(pOwn + pNei);
    }
    else
    {
        phi = this->save(facei, patchi, UNei, Uf_) & normal;
        scalar alphaRhoPhi = alphaNei*rhoNei*UvNei;
        alphaRhoUPhi = alphaRhoPhi*UNei + alphaNei*pNei*normal;
        alphaRhoEPhi = alphaRhoPhi*HNei;

        forAll(alphasOwn, phasei)
        {
            alphaPhis[phasei] = alphasNei[phasei]*UvNei;
            alphaRhoPhis[phasei] = alphaPhis[phasei]*rhosNei[phasei];
        }

        alpha = alphaNei;
        U = UNei;
        p = pNei;
    }

    this->save(facei, patchi, alpha, alphaf_);
    this->save(facei, patchi, U, Uf_);
    this->save(facei, patchi, p, pf_);

    phi *= magSf;
    alphaRhoUPhi *= magSf;
    alphaRhoEPhi *= magSf;
    alphaRhoEPhi += vMesh*magSf*alpha*p;
    forAll(alphasOwn, phasei)
    {
        alphaPhis[phasei] *= magSf;
        alphaRhoPhis[phasei] *= magSf;
    }
}


Foam::scalar Foam::phaseFluxSchemes::HLL::interpolate
(
    const scalar& fOwn, const scalar& fNei,
    const bool rho,
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

// ************************************************************************* //
