/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2019-10-21  Jeff Heylmun:   Moved from rhoCentralFoam to runtime selectable
                            method.
-------------------------------------------------------------------------------License
    This file is part of OpenFOAM.

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

#include "Tadmor.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluxSchemes
{
    defineTypeNameAndDebug(Tadmor, 0);
    addToRunTimeSelectionTable(fluxScheme, Tadmor, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxSchemes::Tadmor::Tadmor
(
    const fvMesh& mesh
)
:
    fluxScheme(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxSchemes::Tadmor::~Tadmor()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluxSchemes::Tadmor::clear()
{
    fluxScheme::clear();
    aPhivOwn_.clear();
    aPhivNei_.clear();
    aSf_.clear();
}


void Foam::fluxSchemes::Tadmor::createSavedFields()
{
    fluxScheme::createSavedFields();
    if (aPhivOwn_.valid())
    {
        return;
    }

    aPhivOwn_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "Tadmor::aPhivOwn",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimVelocity*dimArea, 0.0)
        )
    );
    aPhivNei_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "Tadmor::aPhivNei",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimVelocity*dimArea, 0.0)
        )
    );
    aSf_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "Kurganov::aSf",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimVelocity*dimArea, 0.0)
        )
    );
}


void Foam::fluxSchemes::Tadmor::calculateFluxes
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

    scalar EOwn = eOwn + 0.5*magSqr(UOwn);
    scalar ENei = eNei + 0.5*magSqr(UNei);

    scalar phivOwn(UOwn & Sf);
    scalar phivNei(UNei & Sf);

    scalar cSfOwn(cOwn*magSf);
    scalar cSfNei(cNei*magSf);

    const scalar vMesh(meshPhi(facei, patchi));
    phivOwn -= vMesh;
    phivNei -= vMesh;

    scalar ap(max(max(phivOwn + cSfOwn, phivNei + cSfNei), 0.0));
    scalar am(min(min(phivOwn - cSfOwn, phivNei - cSfNei), 0.0));

    scalar amaxSf(max(mag(am), mag(ap)));
    scalar aSf(-0.5*amaxSf);

    phivOwn *= 0.5;
    phivNei *= 0.5;

    scalar aphivOwn(phivOwn - aSf);
    scalar aphivNei(phivNei + aSf);

    this->save(facei, patchi, aphivOwn, aPhivOwn_);
    this->save(facei, patchi, aphivNei, aPhivNei_);
    this->save(facei, patchi, aSf, aSf_);

    this->save(facei, patchi, 0.5*(UOwn + UNei), Uf_);
    phi = aphivOwn + aphivNei;
    rhoPhi = aphivOwn*rhoOwn + aphivNei*rhoNei;

    rhoUPhi =
    (
        (aphivOwn*rhoOwn*UOwn + aphivNei*rhoNei*UNei)
      + 0.5*(pOwn + pNei)*Sf
    );

    rhoEPhi =
    (
        aphivOwn*(rhoOwn*EOwn + pOwn)
      + aphivNei*(rhoNei*ENei + pNei)
      + aSf*pOwn - aSf*pNei
      + vMesh*0.5*(pOwn + pNei)
    );
}


void Foam::fluxSchemes::Tadmor::calculateFluxes
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

    scalar EOwn = eOwn + 0.5*magSqr(UOwn);
    scalar ENei = eNei + 0.5*magSqr(UNei);

    scalar phivOwn(UOwn & Sf);
    scalar phivNei(UNei & Sf);

    scalar cSfOwn(cOwn*magSf);
    scalar cSfNei(cNei*magSf);

    const scalar vMesh(meshPhi(facei, patchi));
    phivOwn -= vMesh;
    phivNei -= vMesh;

    scalar aOwn(max(max(phivOwn + cSfOwn, phivNei + cSfNei), 0.0));
    scalar aNei(min(min(phivOwn - cSfOwn, phivNei - cSfNei), 0.0));

    scalar amaxSf(max(mag(aNei), mag(aOwn)));
    scalar aSf(-0.5*amaxSf);

    phivOwn *= 0.5;
    phivNei *= 0.5;

    scalar aphivOwn(phivOwn - aSf);
    scalar aphivNei(phivNei + aSf);

    this->save(facei, patchi, aphivOwn, aPhivOwn_);
    this->save(facei, patchi, aphivNei, aPhivNei_);
    this->save(facei, patchi, aSf, aSf_);

    this->save(facei, patchi, 0.5*(UOwn + UNei), Uf_);

    phi = aphivOwn + aphivNei;

    forAll(alphasOwn, phasei)
    {
        alphaPhis[phasei] =
            aphivOwn*alphasOwn[phasei] + aphivNei*alphasNei[phasei];
        alphaRhoPhis[phasei] =
            aphivOwn*alphasOwn[phasei]*rhosOwn[phasei]
          + aphivNei*alphasNei[phasei]*rhosNei[phasei];
    }

    rhoUPhi =
    (
        (aphivOwn*rhoOwn*UOwn + aphivNei*rhoNei*UNei)
      + 0.5*(pOwn + pNei)*Sf
    );

    rhoEPhi =
    (
        aphivOwn*(rhoOwn*EOwn + pOwn)
      + aphivNei*(rhoNei*ENei + pNei)
      + aSf*pOwn - aSf*pNei
      + vMesh*0.5*(pOwn + pNei)
    );
}


Foam::scalar Foam::fluxSchemes::Tadmor::energyFlux
(
    const scalar& rhoOwn, const scalar& rhoNei,
    const vector& UOwn, const vector& UNei,
    const scalar& eOwn, const scalar& eNei,
    const scalar& pOwn, const scalar& pNei,
    const label facei, const label patchi
) const
{
    scalar aphivOwn = getValue(facei, patchi, aPhivOwn_());
    scalar aphivNei = getValue(facei, patchi, aPhivNei_());
    scalar aSf = getValue(facei, patchi, aSf_());

    scalar EOwn = eOwn + 0.5*magSqr(UOwn);
    scalar ENei = eNei + 0.5*magSqr(UNei);

    return
    (
        aphivOwn*(rhoOwn*EOwn + pOwn)
      + aphivNei*(rhoNei*ENei + pNei)
      + aSf*pOwn - aSf*pNei
      + meshPhi(facei, patchi)*0.5*(pOwn + pNei)
    );
}


Foam::scalar Foam::fluxSchemes::Tadmor::interpolate
(
    const scalar& fOwn, const scalar& fNei,
    const label facei, const label patchi
) const
{
    return 0.5*(fOwn + fNei);
}

// ************************************************************************* //
