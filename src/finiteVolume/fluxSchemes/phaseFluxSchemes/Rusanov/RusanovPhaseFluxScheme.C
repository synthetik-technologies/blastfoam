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

#include "RusanovPhaseFluxScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseFluxSchemes
{
    defineTypeNameAndDebug(Rusanov, 0);
    addToRunTimeSelectionTable(phaseFluxScheme, Rusanov, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseFluxSchemes::Rusanov::Rusanov
(
    const surfaceScalarField& phi,
    const scalar residualAlpha
)
:
    phaseFluxScheme(phi, residualAlpha)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseFluxSchemes::Rusanov::~Rusanov()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::phaseFluxSchemes::Rusanov::clear()
{
    phaseFluxScheme::clear();
    lambda_.clear();
}


void Foam::phaseFluxSchemes::Rusanov::createSavedFields()
{
    phaseFluxScheme::createSavedFields();
    if (lambda_.valid())
    {
        lambda_.ref() = Zero;
        return;
    }

    lambda_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                fieldName("lambda"),
                mesh_.time().name(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimVelocity, 0.0)
        )
    );
}


void Foam::phaseFluxSchemes::Rusanov::calculateFluxes
(
    const scalar& alphaOwn, const scalar& alphaNei,
    const scalar& rhoOwn, const scalar& rhoNei,
    const vector& UOwn, const vector& UNei,
    const scalar& eOwn, const scalar& eNei,
    const scalar& pOwn, const scalar& pNei,
    const scalar& cOwn, const scalar& cNei,
    const vector& Sf,
    scalar& phi,
    scalar& alphaRhoPhi,
    vector& alphaRhoUPhi,
    scalar& alphaRhoEPhi,
    const label facei, const label patchi
)
{
    const scalar magSf = mag(Sf);

    const scalar alphaRhoOwn = alphaOwn*rhoOwn;
    const scalar alphaRhoNei = alphaNei*rhoNei;

    const scalar EOwn = eOwn + 0.5*magSqr(UOwn);
    const scalar ENei = eNei + 0.5*magSqr(UNei);

    scalar phivOwn(UOwn & Sf);
    scalar phivNei(UNei & Sf);

    const scalar phiMesh = meshPhi(facei, patchi);
    phivOwn -= phiMesh;
    phivNei -= phiMesh;

    const scalar lambda = max(mag(phivOwn) + cOwn, mag(phivNei) + cNei);
    const scalar phiLambda = lambda*magSf;

    this->save(facei, patchi, lambda, lambda_);
    this->save(facei, patchi, 0.5*(alphaOwn + alphaNei), alphaf_);
    this->save(facei, patchi, 0.5*(UOwn + UNei), Uf_);



    phi = 0.5*(phivOwn + phivNei);

    alphaRhoPhi =
        0.5*(phivOwn*alphaRhoOwn + phivNei*alphaRhoNei)
      - phiLambda*(alphaRhoNei - alphaRhoOwn);

    alphaRhoUPhi =
    (
        0.5
       *(
            phivOwn*alphaOwn*rhoOwn*UOwn
          + phivNei*alphaNei*rhoNei*UNei
          + (alphaOwn*pOwn + alphaNei*pNei)*Sf
        )
      - phiLambda*(alphaRhoNei*UNei - alphaRhoOwn*UOwn)
    );

    alphaRhoEPhi =
    (
        0.5
       *(
            phivOwn*(alphaRhoOwn*EOwn + alphaOwn*pOwn)
          + phivNei*(alphaRhoNei*ENei + alphaNei*pNei)
          + phiMesh*(alphaOwn*pOwn + alphaNei*pNei)
        )
      - phiLambda*(alphaRhoNei*ENei - alphaRhoOwn*EOwn)
    );
}


Foam::scalar Foam::phaseFluxSchemes::Rusanov::calculateAlphaCorrector
(
    const scalar& alphaOwn, const scalar& alphaNei,
    const label facei, const label patchi
) const
{
    return -getValue(facei, patchi, lambda_)*(alphaNei - alphaOwn);
}


Foam::scalar Foam::phaseFluxSchemes::Rusanov::calculateFlux
(
    const scalar& fOwn, const scalar& fNei,
    const scalar& phi,
    const label facei, const label patchi
) const
{
    return
        0.5*phi*(fOwn + fNei)
      - getValue(facei, patchi, lambda_)
       *(fNei - fOwn)
       *getValue(facei, patchi, mesh_.magSf());
}


Foam::scalar Foam::phaseFluxSchemes::Rusanov::interpolate
(
    const scalar& fOwn, const scalar& fNei,
    const label facei, const label patchi
) const
{
    return 0.5*(fOwn + fNei);
}

// ************************************************************************* //
