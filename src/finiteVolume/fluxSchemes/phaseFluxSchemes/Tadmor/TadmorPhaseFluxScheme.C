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

#include "TadmorPhaseFluxScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseFluxSchemes
{
    defineTypeNameAndDebug(Tadmor, 0);
    addToRunTimeSelectionTable(phaseFluxScheme, Tadmor, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseFluxSchemes::Tadmor::Tadmor
(
    const surfaceScalarField& phi
)
:
    phaseFluxScheme(phi)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseFluxSchemes::Tadmor::~Tadmor()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::phaseFluxSchemes::Tadmor::clear()
{
    phaseFluxScheme::clear();
    aPhivOwn_.clear();
    aPhivNei_.clear();
}


void Foam::phaseFluxSchemes::Tadmor::createSavedFields()
{
    phaseFluxScheme::createSavedFields();

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
}


void Foam::phaseFluxSchemes::Tadmor::calculateFluxes
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

    this->save(facei, patchi, 0.5*(UOwn + UNei), Uf_);

    phi = aphivOwn + aphivNei;

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


Foam::scalar Foam::phaseFluxSchemes::Tadmor::calculateFlux
(
    const scalar& fOwn, const scalar& fNei,
    const scalar& phi,
    const label facei, const label patchi
) const
{
    return
        getValue(facei, patchi, aPhivOwn_)*fOwn
      + getValue(facei, patchi, aPhivNei_)*fNei;
}


Foam::scalar Foam::phaseFluxSchemes::Tadmor::interpolate
(
    const scalar& fOwn, const scalar& fNei,
    const bool rho,
    const label facei, const label patchi
) const
{
    const scalar aphivOwn(getValue(facei, patchi, aPhivOwn_));
    const scalar aphivNei(getValue(facei, patchi, aPhivNei_));
    const scalar aphiv(aphivOwn + aphivNei);
    if (mag(aphiv) > small)
    {
        return (fOwn*aphivOwn + fNei*aphivNei)/aphiv;
    }
    return 0.5*(fOwn + fNei);
}

// ************************************************************************* //
