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
    const scalar& rhoOwn, const scalar& rhoNei,
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
    scalar HOwn(EOwn + pOwn/rhoOwn);

    scalar ENei = eNei + 0.5*magSqr(UNei);
    scalar HNei(ENei + pNei/rhoNei);

    const scalar vMesh(meshPhi(facei, patchi)/magSf);
    scalar UvOwn((UOwn & normal) - vMesh);
    scalar UvNei((UNei & normal) - vMesh);

    scalar SOwn(stabilise(min(UvOwn - cOwn, UvNei - cNei), small));
    scalar SNei(stabilise(max(UvOwn + cOwn, UvNei + cNei), small));

    this->save(facei, patchi, SOwn, SOwn_);
    this->save(facei, patchi, SNei, SNei_);
    this->save(facei, patchi, UvOwn, UvOwn_);
    this->save(facei, patchi, UvNei, UvNei_);

    scalar alphaP;
    if (SOwn >= 0)
    {
        phi = this->save(facei, patchi, UOwn, Uf_) & normal;
        alphaPhi = alphaOwn*UvOwn;
        alphaRhoPhi = alphaOwn*rhoOwn*UvOwn;
        alphaRhoUPhi = alphaOwn*UOwn*rhoOwn*UvOwn + alphaOwn*pOwn*normal;
        alphaRhoEPhi = alphaOwn*UvOwn*rhoOwn*HOwn;
        alphaP = alphaOwn*pOwn;
    }
    else if (SOwn < 0 && SNei >= 0)
    {
        scalar alphaRhoOwn = alphaOwn*rhoOwn;
        scalar alphaRhoNei = alphaNei*rhoNei;
        vector alphaRhoUOwn = alphaRhoOwn*UOwn;
        vector alphaRhoUNei = alphaRhoNei*UNei;
        scalar alphaRhoPhiOwn = alphaRhoOwn*UvOwn;
        scalar alphaRhoPhiNei = alphaRhoNei*UvNei;
        vector alphaRhoUPhiOwn = alphaRhoUOwn*UvOwn + pOwn*normal;
        vector alphaRhoUPhiNei = alphaRhoUNei*UvNei + pNei*normal;
        phi =
            this->save
            (
                facei,
                patchi,
                (
                    SNei*alphaRhoUNei - SOwn*alphaRhoUOwn
                  + alphaRhoUPhiOwn - alphaRhoUPhiNei
                )
               /(
                   SNei*alphaRhoNei - SOwn*alphaRhoOwn
                 + alphaRhoPhiNei - alphaRhoPhiOwn
                ),
                Uf_
            ) & normal;

        alphaPhi =
            (
                SNei*alphaOwn - SOwn*alphaNei
              + SOwn*SNei*(alphaNei - alphaOwn)
            )/(SNei - SOwn);

        alphaRhoPhi =
            (
                SNei*alphaRhoPhiOwn - SOwn*alphaRhoPhiNei
              + SOwn*SNei*(alphaRhoNei - alphaRhoOwn)
            )/(SNei - SOwn);

        alphaRhoUPhi =
            (
                SNei*alphaRhoUPhiOwn - SOwn*alphaRhoUPhiNei
              + SOwn*SNei*(alphaRhoUNei - alphaRhoUOwn)
            )/(SNei - SOwn);

        alphaRhoEPhi =
            (
                SNei*alphaRhoOwn*HOwn*UvOwn
              - SOwn*alphaRhoNei*HNei*UvNei
              + SOwn*SNei*(alphaRhoNei*ENei - alphaRhoOwn*EOwn)
            )/(SNei - SOwn);
        alphaP = (SNei*alphaOwn*pOwn - SOwn*alphaNei*pNei)/(SNei - SOwn);
    }
    else
    {
        phi = this->save(facei, patchi, UNei, Uf_) & normal;
        alphaPhi = alphaNei*UvNei;
        alphaRhoPhi = alphaNei*rhoNei*UvNei;
        alphaRhoUPhi = alphaNei*UNei*rhoNei*UvNei + alphaNei*pNei*normal;
        alphaRhoEPhi = alphaNei*rhoNei*HNei*UvNei;
        alphaP = alphaNei*pNei;
    }
    phi *= magSf;
    alphaPhi *= magSf;
    alphaRhoPhi *= magSf;
    alphaRhoUPhi *= magSf;
    alphaRhoEPhi *= magSf;
    alphaRhoEPhi += vMesh*magSf*alphaP;
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
    NotImplemented;
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
