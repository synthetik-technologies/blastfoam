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

#include "KurganovPhaseFluxScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseFluxSchemes
{
    defineTypeNameAndDebug(Kurganov, 0);
    addToRunTimeSelectionTable(phaseFluxScheme, Kurganov, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseFluxSchemes::Kurganov::Kurganov(const surfaceScalarField& phi)
:
    Tadmor(phi)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseFluxSchemes::Kurganov::~Kurganov()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::phaseFluxSchemes::Kurganov::clear()
{
    Tadmor::clear();
    aOwn_.clear();
    aNei_.clear();
}


void Foam::phaseFluxSchemes::Kurganov::createSavedFields()
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
                mesh_.time().name(),
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
                mesh_.time().name(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimless, 0.0)
        )
    );
}

void Foam::phaseFluxSchemes::Kurganov::calculateFluxes
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

    scalar EOwn = eOwn + 0.5*magSqr(UOwn);
    scalar ENei = eNei + 0.5*magSqr(UNei);

    scalar phivOwn(UOwn & Sf);
    scalar phivNei(UNei & Sf);

    scalar cSfOwn(cOwn*magSf);
    scalar cSfNei(cNei*magSf);

    const scalar vMesh(meshPhi(facei, patchi));
    phivOwn -= vMesh;
    phivNei -= vMesh;

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

    this->save(facei, patchi, aOwn*UOwn + aNei*UNei, Uf_);

    alphaPhi = aphivOwn*alphaOwn + aphivNei*alphaNei;
    alphaRhoPhi = aphivOwn*alphaOwn*rhoOwn + aphivNei*alphaNei*rhoNei;

    alphaRhoUPhi =
    (
        (
            aphivOwn*alphaOwn*rhoOwn*UOwn
          + aphivNei*alphaNei*rhoNei*UNei
        )
      + 0.5*(alphaOwn*pOwn + alphaNei*pNei)*Sf
    );

    alphaRhoEPhi =
    (
        aphivOwn*(alphaOwn*(rhoOwn*EOwn + pOwn))
      + aphivNei*(alphaNei*(rhoNei*ENei + pNei))
      + aSf*(alphaOwn*pOwn - alphaNei*pNei)
      + vMesh*0.5*(alphaOwn*pOwn + alphaNei*pNei)
    );
}


Foam::scalar Foam::phaseFluxSchemes::Kurganov::interpolate
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
