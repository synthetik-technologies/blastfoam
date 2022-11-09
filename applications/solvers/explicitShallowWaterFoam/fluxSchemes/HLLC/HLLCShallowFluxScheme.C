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

#include "HLLCShallowFluxScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace shallowFluxSchemes
{
    defineTypeNameAndDebug(HLLC, 0);
    addToRunTimeSelectionTable(shallowFluxScheme, HLLC, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shallowFluxSchemes::HLLC::HLLC
(
    surfaceScalarField& phi,
    surfaceScalarField& hPhi,
    surfaceVectorField& hUPhi,
    const dimensionedVector& g
)
:
    shallowFluxScheme(phi, hPhi, hUPhi, g)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::shallowFluxSchemes::HLLC::~HLLC()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::shallowFluxSchemes::HLLC::clear()
{
    shallowFluxScheme::clear();
}


void Foam::shallowFluxSchemes::HLLC::createSavedFields()
{
    shallowFluxScheme::createSavedFields();
}


void Foam::shallowFluxSchemes::HLLC::calculateFluxes
(
    const scalar& hOwn, const scalar& hNei,
    const scalar& h0Own, const scalar& h0Nei,
    const vector& UOwn, const vector& UNei,
    const vector& Sf,
    scalar& phi,
    scalar& hPhi,
    vector& hUPhi,
    const label facei, const label patchi
)
{
    scalar magSf = mag(Sf);
    vector normal = Sf/magSf;

    const scalar vMesh(meshPhi(facei, patchi)/magSf);

    scalar UvOwn = (UOwn & normal) - vMesh;
    scalar UvNei = (UNei & normal) - vMesh;

    scalar cOwn = wavespeed(hOwn);
    scalar cNei = wavespeed(hNei);

    scalar  hStar =
        0.5*(hOwn + hNei)
      - 0.25*(UvNei - UvOwn)*(hOwn + hNei)/(cOwn + cNei);

    scalar qOwn  =
        hStar > hOwn
      ? sqrt(0.5*((hStar + hOwn)*hStar)/sqr(max(hOwn, 1e-10)))
      : 1.0;
    scalar qNei  =
        hStar > hNei
      ? sqrt(0.5*((hStar + hNei*hStar)/sqr(max(hNei, 1e-10))))
      : 1.0;

    scalar SOwn(UvOwn - cOwn*qOwn);
    scalar SNei(UvNei + cNei*qNei);
    scalar  SStar =
        0.5*(UvOwn + UvNei)
      - 0.25*(hNei - hOwn)*(cOwn + cNei)/(hOwn + hNei);

    scalar hHLLC =
        (SNei*hNei - SOwn*hOwn - (hNei*UvNei - hOwn*UvOwn))
       /(SNei - SOwn);

    scalar qHLLC =
        (
            SNei*hNei*UvNei
          - SOwn*hOwn*UvOwn
          - (hNei*sqr(UvNei) + 0.5*magg_*sqr(hNei))
          + (hOwn*sqr(UvOwn) + 0.5*magg_*sqr(hOwn))
        )/(SNei - SOwn);

    scalar dh0 = h0Nei - h0Own;
    if (SOwn >= 0)
    {
        phi = UvOwn;
        hPhi = phi*hOwn;
        hUPhi = hPhi*UOwn + 0.5*magg_*sqr(hOwn)*normal;
    }
    else if (SOwn < 0 && SStar >= 0)
    {
        scalar f = (SOwn - UvOwn)/(SOwn - SStar);

        phi = f*SStar;

        hPhi = hOwn*phi;

        vector hUStar = f*hOwn*(UOwn - (UvOwn - SStar)*normal);
        hUPhi =
            hOwn*UvOwn*UOwn + 0.5*magg_*sqr(hOwn)*normal
          + SOwn*(hUStar - hOwn*UOwn);
    }
    else if (SStar < 0 && SNei >= 0)
    {
        scalar f = (SNei - UvNei)/(SNei - SStar);

        phi = f*SStar;

        hPhi = hNei*phi;

        vector hUStar = f*hNei*(UNei - (UvNei - SStar)*normal);
        hUPhi =
            hNei*UvNei*UNei + 0.5*magg_*sqr(hNei)*normal
          + SNei*(hUStar - hNei*UNei);
    }
    else
    {
        phi = UvNei;
        hPhi = phi*hNei;
        hUPhi = hPhi*UNei + 0.5*magg_*sqr(hNei)*normal;
    }
    phi *= magSf;
    hPhi *= magSf;
    hUPhi *= magSf;
}

// ************************************************************************* //
