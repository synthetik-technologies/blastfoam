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

#include "RoeFluxScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluxSchemes
{
    defineTypeNameAndDebug(Roe, 0);
    addToRunTimeSelectionTable(fluxScheme, Roe, singlePhase);
//     addToRunTimeSelectionTable(fluxScheme, Roe, multiphase);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxSchemes::Roe::Roe
(
    const fvMesh& mesh
)
:
    fluxScheme(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxSchemes::Roe::~Roe()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluxSchemes::Roe::clear()
{
    fluxScheme::clear();
}

void Foam::fluxSchemes::Roe::createSavedFields()
{
    fluxScheme::createSavedFields();
}


void Foam::fluxSchemes::Roe::calculateFluxes
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

    scalar HOwn = eOwn + 0.5*magSqr(UOwn) + pOwn/rhoOwn;
    scalar HNei = eNei + 0.5*magSqr(UNei) + pNei/rhoNei;

    const scalar vMesh(meshPhi(facei, patchi)/magSf);
    scalar UvOwn((UOwn & normal) - vMesh);
    scalar UvNei((UNei & normal) - vMesh);

    scalar wOwn(sqrt(rhoOwn)/(sqrt(rhoOwn) + sqrt(rhoNei)));
    scalar wNei(1.0 - wOwn);

    scalar rhoTilde(sqrt(rhoOwn*rhoNei));

    vector UTilde(UOwn*wOwn + UNei*wNei);
    scalar UvTilde(UTilde & normal);

    scalar HTilde(HOwn*wOwn + HNei*wNei);

    scalar cTilde(cOwn*wOwn + cNei*wNei);

    scalar deltaRho(rhoNei - rhoOwn);
    vector deltaU(UNei - UOwn);
    scalar deltaUv(deltaU & normal);
    scalar deltaP(pNei - pOwn);

    scalar lambda1(mag(UvTilde));
    scalar lambda2(mag(UvTilde + cTilde));
    scalar lambda3(mag(UvTilde - cTilde));

    scalar alpha1(deltaRho - deltaP/sqr(cTilde));
    scalar alpha2((deltaP + rhoTilde*cTilde*deltaUv)/(2.0*sqr(cTilde)));
    scalar alpha3((deltaP - rhoTilde*cTilde*deltaUv)/(2.0*sqr(cTilde)));

    // U Row
    vector K21(UTilde);
    vector K224((UTilde + cTilde*normal));
    vector K25((UTilde - cTilde*normal));

    // E row
    scalar K31(0.5*magSqr(UTilde));
    scalar K324((HTilde + cTilde*UvTilde));
    scalar K35((HTilde - cTilde*UvTilde));

//     {
//         scalar eps = 0.1*cTilde; //adjustable parameter
//
//         if (lambda1 < eps || lambda2 < eps || lambda3 < eps)
//         {
//             lambda1 = (sqr(lambda1) + sqr(eps))/(2.0*eps);
//             lambda2 = (sqr(lambda2) + sqr(eps))/(2.0*eps);
//             lambda3 = (sqr(lambda3) + sqr(eps))/(2.0*eps);
//         }
//
//         // First eigenvalue: U - c
//         eps = 2.0*max(0.0, (UvNei - cNei) - (UvOwn - cOwn));
//         if (lambda1 < eps)
//         {
//             lambda1 = (sqr(lambda1) + sqr(eps))/(2.0*eps);
//         }
//
//         // Second eigenvalue: U
//         eps = 2.0*max(0.0, UvNei - UvOwn);
//         if (lambda2 < eps)
//         {
//             lambda2 = (sqr(lambda2) + sqr(eps))/(2.0*eps);
//         }
//
//         // Third eigenvalue: U + c
//         eps = 2.0*max(0.0, (UvNei + cNei) - (UvOwn + cOwn));
//         if (lambda3 < eps)
//         {
//             lambda3 = (sqr(lambda3) + sqr(eps))/(2.0*eps);
//         }
//     }

//     this->save(facei, patchi, lambda1, lambda1_);
//     this->save(facei, patchi, lambda2, lambda2_);
//     this->save(facei, patchi, lambda3, lambda3_);


    scalar rhoPhiOwn(rhoOwn*UvOwn);
    scalar rhoPhiNei(rhoNei*UvNei);

    vector rhoUPhiOwn(UOwn*rhoPhiOwn + pOwn*normal);
    vector rhoUPhiNei(UNei*rhoPhiNei + pNei*normal);

    scalar rhoEPhiOwn(HOwn*rhoPhiOwn);
    scalar rhoEPhiNei(HNei*rhoPhiNei);

    // Compute fluxes
    rhoPhi =
        0.5*magSf
       *(
            rhoPhiOwn + rhoPhiNei
          - (
                lambda1*alpha1
              + lambda2*alpha2
              + lambda3*alpha3
            )
       );
    phi = rhoPhi/rhoTilde;

    rhoUPhi =
        0.5*magSf
       *(
            rhoUPhiOwn + rhoUPhiNei
          - (
                lambda1*alpha1*K21
              + lambda2*alpha2*K224
              + lambda3*alpha3*K25
            )
       );

    rhoEPhi =
        0.5*magSf
       *(
            rhoEPhiOwn + rhoEPhiNei
          - (
                lambda1*alpha1*K31
              + lambda2*alpha2*K324
              + lambda3*alpha3*K35
            )
       );
}


void Foam::fluxSchemes::Roe::calculateFluxes
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
   NotImplemented;
}


Foam::scalar Foam::fluxSchemes::Roe::energyFlux
(
    const scalar& rhoOwn, const scalar& rhoNei,
    const vector& UOwn, const vector& UNei,
    const scalar& eOwn, const scalar& eNei,
    const scalar& pOwn, const scalar& pNei,
    const vector& Sf,
    const label facei, const label patchi
) const
{
    NotImplemented;
    return 0.0;
}


Foam::scalar Foam::fluxSchemes::Roe::interpolate
(
    const scalar& fOwn, const scalar& fNei,
    const bool isDensity,
    const label facei, const label patchi
) const
{
    NotImplemented;
    return 0.0;
}


// ************************************************************************* //
