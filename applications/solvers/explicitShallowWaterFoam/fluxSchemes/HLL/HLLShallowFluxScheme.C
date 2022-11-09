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

#include "HLLShallowFluxScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace shallowFluxSchemes
{
    defineTypeNameAndDebug(HLL, 0);
    addToRunTimeSelectionTable(shallowFluxScheme, HLL, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shallowFluxSchemes::HLL::HLL
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

Foam::shallowFluxSchemes::HLL::~HLL()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::shallowFluxSchemes::HLL::clear()
{
    shallowFluxScheme::clear();
}


void Foam::shallowFluxSchemes::HLL::createSavedFields()
{
    shallowFluxScheme::createSavedFields();
}


void Foam::shallowFluxSchemes::HLL::calculateFluxes
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

    scalar SOwn(min(UvOwn - cOwn, UvNei - cNei));
    scalar SNei(max(UvOwn + cOwn, UvNei + cNei));

    if (SOwn >= 0)
    {
        phi = magSf*UvOwn;
        hPhi = phi*hOwn;
        hUPhi = hPhi*UOwn + 0.5*magg_*sqr(hOwn)*Sf;
    }
    else if (SOwn < 0 && SNei >= 0)
    {
        vector hUOwn = hOwn*UOwn;
        vector hUNei = hNei*UNei;
        scalar hPhiOwn = hOwn*UvOwn;
        scalar hPhiNei = hNei*UvNei;
        vector hUPhiOwn = hPhiOwn*UOwn + 0.5*magg_*sqr(hOwn)*normal;
        vector hUPhiNei = hPhiNei*UNei + 0.5*magg_*sqr(hNei)*normal;

        // hU_HLL / h_HLL
        phi =
            (
                (SNei*hUNei - SOwn*hUOwn + hUPhiOwn - hUPhiNei)
               /(SNei*hNei - SOwn*hOwn + hPhiOwn - hPhiNei)
            ) & Sf;

        hPhi =
            (
                SNei*hPhiOwn - SOwn*hPhiNei
              + SOwn*SNei*(hNei - hOwn)
            )/(SNei - SOwn)*magSf;

        hUPhi =
            (
                SNei*hUPhiOwn - SOwn*hUPhiNei
              + SOwn*SNei*(hUNei - hUOwn)
            )/(SNei - SOwn)*magSf;
    }
    else
    {
        phi = magSf*UvNei;
        hPhi = phi*hNei;
        hUPhi = hPhi*UNei + 0.5*magg_*sqr(hNei)*Sf;
    }
}

// ************************************************************************* //
