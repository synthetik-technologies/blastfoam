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

#include "HLL.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluxSchemes
{
    defineTypeNameAndDebug(HLL, 0);
    addToRunTimeSelectionTable(fluxScheme, HLL, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxSchemes::HLL::HLL
(
    const fvMesh& mesh
)
:
    fluxScheme(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxSchemes::HLL::~HLL()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluxSchemes::HLL::calculateFluxes
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

    scalar wOwn(sqrt(rhoOwn)/(sqrt(rhoOwn) + sqrt(rhoNei)));
    scalar wNei(1.0 - wOwn);

    scalar cTilde(cOwn*wOwn + cNei*wNei);
    scalar UvTilde(UvOwn*wOwn + UvNei*wNei);

    scalar SOwn(min(UvOwn - cOwn, UvTilde - cTilde));
    scalar SNei(max(UvNei + cNei, UvTilde + cTilde));

    if (SOwn >= 0)
    {
        phi = UvOwn;
        rhoUPhi = UOwn*rhoOwn*UvOwn + pOwn*normal;
        rhoEPhi = UvOwn*rhoOwn*HOwn;

        forAll(alphasOwn, phasei)
        {
            alphaPhis[phasei] = alphasOwn[phasei]*UvOwn;
            alphaRhoPhis[phasei] = alphasOwn[phasei]*rhosOwn[phasei]*UvOwn;
        }
    }
    else if (SOwn < 0 && SNei >= 0)
    {
        phi = (SNei*UvOwn - SOwn*UvNei)/(SNei - SOwn);

        vector rhoUOwn = rhoOwn*UOwn;
        vector rhoUNei = rhoNei*UNei;
        rhoUPhi =
            (
                SNei*(rhoUOwn*UvOwn + pOwn*normal)
              - SOwn*(rhoUNei*UvNei + pNei*normal)
              + SOwn*SNei*(rhoUNei - rhoUOwn)
            )/(SNei - SOwn);

        rhoEPhi =
            (
                SNei*rhoOwn*HOwn*UvOwn
              - SOwn*rhoNei*HNei*UvNei
              + SOwn*SNei
               *(
                    rhoNei*ENei
                  - rhoOwn*EOwn
                )
            )/(SNei - SOwn);

        forAll(alphasOwn, phasei)
        {
            alphaPhis[phasei] =
            (
                SNei*alphasOwn[phasei]*UvOwn - SOwn*alphasNei[phasei]*UvNei
              + SOwn*SNei*(alphasNei[phasei] - alphasOwn[phasei])
            )/(SNei - SOwn);

            scalar alphaRhoOwn = alphasOwn[phasei]*rhosOwn[phasei];
            scalar alphaRhoNei = alphasNei[phasei]*rhosNei[phasei];
            alphaRhoPhis[phasei] =
                (
                    SNei*alphaRhoOwn*UvOwn - SOwn*alphaRhoNei*UvNei
                  + SOwn*SNei*(alphaRhoNei - alphaRhoOwn)
                )/(SNei - SOwn);
        }
    }
    else if (SNei < 0)
    {
        phi = UvNei;
        rhoUPhi = UNei*rhoNei*UvNei + pNei*normal;
        rhoEPhi = rhoNei*HNei*UvNei;

        forAll(alphasOwn, phasei)
        {
            alphaPhis[phasei] = alphasNei[phasei]*UvNei;
            alphaRhoPhis[phasei] = alphasNei[phasei]*rhosNei[phasei]*UvNei;
        }
    }
    phi *= magSf;
    rhoUPhi *= magSf;
    rhoEPhi *= magSf;
    forAll(alphasOwn, phasei)
    {
        alphaPhis[phasei] *= magSf;
        alphaRhoPhis[phasei] *= magSf;
    }
}

// ************************************************************************* //
