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

#include "AUSMPlus.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluxSchemes
{
    defineTypeNameAndDebug(AUSMPlus, 0);
    addToRunTimeSelectionTable(fluxScheme, AUSMPlus, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxSchemes::AUSMPlus::AUSMPlus(const fvMesh& mesh)
:
    fluxScheme(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxSchemes::AUSMPlus::~AUSMPlus()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fluxSchemes::AUSMPlus::calculateFluxes
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
    scalar& rhoEPhi
)
{
    scalar magSf = mag(Sf);
    vector normal = Sf/magSf;

    scalar EOwn = eOwn + 0.5*magSqr(UOwn);
    scalar HOwn(EOwn + pOwn/rhoOwn);

    scalar ENei = eNei + 0.5*magSqr(UNei);
    scalar HNei(ENei + pNei/rhoNei);

    scalar UvOwn(UOwn & normal);
    scalar UvNei(UNei & normal);

    scalar c12(0.5*(cOwn + cNei));

    // Compute slpit Mach numbers
    scalar MaOwn(UvOwn/c12);
    scalar MaNei(UvNei/c12);
    scalar magMaOwn(mag(MaOwn));
    scalar magMaNei(mag(MaNei));

    scalar Ma4Own = max(MaOwn, 0.0);
    if (magMaOwn < 1)
    {
        Ma4Own = 0.25*sqr(MaOwn + 1.0) + beta_*sqr(sqr(MaOwn) - 1.0);
    }

    scalar Ma4Nei = min(MaNei, 0.0);
    if (magMaNei < 1)
    {
        Ma4Nei = -0.25*sqr(MaNei - 1.0) - beta_*sqr(sqr(MaNei) - 1.0);
    }
    scalar Ma12(Ma4Own + Ma4Nei);

    scalar P5Own
    (
        pos0(magMaOwn - 1.0)*pos(MaOwn)
      + neg(magMaOwn - 1.0)
       *(
            0.25*sqr(MaOwn + 1.0)*(2.0 - MaOwn)
          + alpha_*MaOwn*sqr(sqr(MaOwn) - 1.0)
        )
    );

    scalar P5Nei
    (
        pos0(magMaNei - 1.0)*neg(MaNei)
      + neg(magMaNei - 1.0)
       *(
            0.25*sqr(MaNei - 1.0)*(2.0 + MaNei)
          - alpha_*MaNei*sqr(sqr(MaNei) - 1.0)
        )
    );

    scalar P12(P5Own*pOwn + P5Nei*pNei);

    phi = magSf*c12*Ma12;

    if (Ma12 >= 0)
    {
        rhoUPhi = rhoOwn*UOwn;
        rhoEPhi = rhoOwn*HOwn;
        forAll(alphasOwn, phasei)
        {
            alphaPhis[phasei] = alphasOwn[phasei];
            alphaRhoPhis[phasei] = alphasOwn[phasei]*rhosOwn[phasei];
        }
    }
    else
    {
        rhoUPhi = rhoNei*UNei;
        rhoEPhi = rhoNei*HNei;
        forAll(alphasOwn, phasei)
        {
            alphaPhis[phasei] = alphasNei[phasei];
            alphaRhoPhis[phasei] = alphasNei[phasei]*rhosNei[phasei];
        }
    }
    rhoUPhi *= phi;
    rhoUPhi += P12*Sf;
    rhoEPhi *= phi;

    forAll(alphasOwn, phasei)
    {
        alphaPhis[phasei] *= phi;
        alphaRhoPhis[phasei] *= phi;
    }
}
