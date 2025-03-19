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

#include "KurganovFluxScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluxSchemes
{
    defineTypeNameAndDebug(Kurganov, 0);
    addToRunTimeSelectionTable(fluxScheme, Kurganov, singlePhase);
    addToRunTimeSelectionTable(fluxScheme, Kurganov, multiphase);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxSchemes::Kurganov::Kurganov(const surfaceScalarField& phi)
:
    Tadmor(phi)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxSchemes::Kurganov::~Kurganov()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluxSchemes::Kurganov::clear()
{
    Tadmor::clear();
    aOwn_.clear();
    aNei_.clear();
}


void Foam::fluxSchemes::Kurganov::createSavedFields()
{
    Tadmor::createSavedFields();
    if (aOwn_.valid())
    {
        return;
    }
    aOwn_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                fieldName("aOwn"),
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimless, 0.0)
        )
    );
    aNei_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                fieldName("aNei"),
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimless, 0.0)
        )
    );
}


void Foam::fluxSchemes::Kurganov::calculateFluxes
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

    const scalar meshFlux = meshPhi(facei, patchi);
    phivOwn -= meshFlux;
    phivNei -= meshFlux;

    scalar ap
    (
        max(max(phivOwn + cSfOwn, phivNei + cSfNei), 0.0)
    );
    scalar am
    (
        min(min(phivOwn - cSfOwn, phivNei - cSfNei), 0.0)
    );

    scalar aOwn(ap/(ap - am));
    scalar aSf(am*aOwn);
    scalar aNei(1.0 - aOwn);

    phivOwn *= aOwn;
    phivNei *= aNei;

    scalar aphivOwn(phivOwn - aSf);
    scalar aphivNei(phivNei + aSf);

    this->save(facei, patchi, aphivOwn, aPhivOwn_);
    this->save(facei, patchi, aphivNei, aPhivNei_);
    this->save(facei, patchi, aOwn, aOwn_);
    this->save(facei, patchi, aNei, aNei_);

    if (needEnergyFlux)
    {
        this->save(facei, patchi, aSf, aSf_);
    }

    phi = aphivOwn + aphivNei;
    this->save
    (
        facei,
        patchi,
        (aOwn*UOwn + aNei*UNei),
        Uf_
    );
    rhoPhi = aphivOwn*rhoOwn + aphivNei*rhoNei;

    rhoUPhi =
    (
        (aphivOwn*rhoOwn*UOwn + aphivNei*rhoNei*UNei)
      + (aOwn*pOwn + aNei*pNei)*Sf
    );

    rhoEPhi =
    (
        aphivOwn*(rhoOwn*EOwn + pOwn)
      + aphivNei*(rhoNei*ENei + pNei)
      + aSf*pOwn - aSf*pNei
      + meshFlux*(aOwn*pOwn + aNei*pNei)
    );
}


Foam::scalar Foam::fluxSchemes::Kurganov::energyFlux
(
    const scalar& rhoOwn, const scalar& rhoNei,
    const vector& UOwn, const vector& UNei,
    const scalar& eOwn, const scalar& eNei,
    const scalar& pOwn, const scalar& pNei,
    const vector& Sf,
    const label facei, const label patchi
) const
{
    scalar aphivOwn = getValue(facei, patchi, aPhivOwn_);
    scalar aphivNei = getValue(facei, patchi, aPhivNei_);
    scalar aOwn = getValue(facei, patchi, aOwn_);
    scalar aNei = getValue(facei, patchi, aNei_);
    scalar aSf = getValue(facei, patchi, aSf_);

    scalar EOwn = eOwn + 0.5*magSqr(UOwn);
    scalar ENei = eNei + 0.5*magSqr(UNei);


    return
    (
        aphivOwn*(rhoOwn*EOwn + pOwn)
      + aphivNei*(rhoNei*ENei + pNei)
      + aSf*pOwn - aSf*pNei
      + meshPhi(facei, patchi)*(aOwn*pOwn + aNei*pNei)
    );
}


Foam::scalar Foam::fluxSchemes::Kurganov::interpolate
(
    const scalar& fOwn, const scalar& fNei,
    const label facei, const label patchi
) const
{
    return
        getValue(facei, patchi, aOwn_)*fOwn
      + getValue(facei, patchi, aNei_)*fNei;
}

// ************************************************************************* //
