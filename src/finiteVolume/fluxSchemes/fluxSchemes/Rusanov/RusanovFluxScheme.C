/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "RusanovFluxScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluxSchemes
{
    defineTypeNameAndDebug(Rusanov, 0);
    addToRunTimeSelectionTable(fluxScheme, Rusanov, singlePhase);
//     addToRunTimeSelectionTable(fluxScheme, Rusanov, multiphase);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxSchemes::Rusanov::Rusanov(const surfaceScalarField& phi)
:
    fluxScheme(phi)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxSchemes::Rusanov::~Rusanov()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluxSchemes::Rusanov::clear()
{
    fluxScheme::clear();
    phivOwn_.clear();
    phivNei_.clear();
    lambda_.clear();
}

void Foam::fluxSchemes::Rusanov::createSavedFields()
{
    fluxScheme::createSavedFields();
    if (phivOwn_.valid())
    {
        return;
    }

    phivOwn_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                fieldName("phivOwn"),
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimVelocity*dimArea, 0.0)
        )
    );
    phivNei_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                fieldName("phivNei"),
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimVelocity*dimArea, 0.0)
        )
    );
    lambda_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                fieldName("lambda"),
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimVelocity*dimArea, 0.0)
        )
    );
}


void Foam::fluxSchemes::Rusanov::calculateFluxes
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
    vector normal = Sf/magSf;

    const scalar vMesh(meshPhi(facei, patchi)/magSf);
    scalar UvOwn((UOwn & normal) - vMesh);
    scalar UvNei((UNei & normal) - vMesh);

    scalar lambda = max(mag(UvOwn) + cOwn, mag(UvNei) + cNei)*magSf;

    this->save(facei, patchi, UvOwn*magSf, phivOwn_);
    this->save(facei, patchi, UvNei*magSf, phivNei_);
    this->save(facei, patchi, lambda, lambda_);

    phi = 0.5*(UvOwn + UvNei);

    scalar rhoPhiOwn = rhoOwn*UvOwn*magSf;
    scalar rhoPhiNei = rhoNei*UvNei*magSf;

    rhoPhi = 0.5*(rhoPhiOwn + rhoPhiNei - lambda*(rhoNei - rhoOwn));

    rhoUPhi =
        0.5
       *(
            rhoPhiOwn*UOwn + rhoPhiNei*UNei
          + (pOwn + pNei)*Sf
          - lambda*(rhoNei*UNei - rhoOwn*UOwn)
        );

    scalar EOwn = eOwn + 0.5*magSqr(UOwn);
    scalar ENei = eNei + 0.5*magSqr(UNei);
    rhoEPhi =
        0.5
       *(
            rhoPhiOwn*(EOwn + pOwn/rhoOwn)
          + rhoPhiNei*(ENei + pNei/rhoNei)
          - lambda*(rhoNei*ENei - rhoOwn*EOwn)
        );
}


Foam::scalar Foam::fluxSchemes::Rusanov::energyFlux
(
    const scalar& rhoOwn, const scalar& rhoNei,
    const vector& UOwn, const vector& UNei,
    const scalar& eOwn, const scalar& eNei,
    const scalar& pOwn, const scalar& pNei,
    const vector& Sf,
    const label facei, const label patchi
) const
{
    scalar phivOwn(getValue(facei, patchi, phivOwn_));
    scalar phivNei(getValue(facei, patchi, phivNei_));

    scalar lambda(getValue(facei, patchi, lambda_));

    scalar EOwn = eOwn + 0.5*magSqr(UOwn);
    scalar ENei = eNei + 0.5*magSqr(UNei);
    return
        0.5
       *(
            phivOwn*rhoOwn*(EOwn + pOwn/rhoOwn)
          + phivNei*rhoNei*(ENei + pNei/rhoNei)
          - lambda*(rhoNei*ENei - rhoOwn*EOwn)
        );
}


Foam::scalar Foam::fluxSchemes::Rusanov::interpolate
(
    const scalar& fOwn, const scalar& fNei,
    const label facei, const label patchi
) const
{
    NotImplemented;
    return 0.0;
}


Foam::scalar Foam::fluxSchemes::Rusanov::calculateFlux
(
    const scalar& fOwn, const scalar& fNei,
    const scalar& phi,
    const label facei, const label patchi
) const
{
    return
        0.5
       *(
            getValue(facei, patchi, phivOwn_)*fOwn
          + getValue(facei, patchi, phivNei_)*fNei
          - getValue(facei, patchi, lambda_)*(fNei - fOwn)
        );
}


// ************************************************************************* //
