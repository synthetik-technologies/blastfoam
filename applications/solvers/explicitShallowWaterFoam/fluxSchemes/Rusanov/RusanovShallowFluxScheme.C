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

#include "RusanovShallowFluxScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace shallowFluxSchemes
{
    defineTypeNameAndDebug(Rusanov, 0);
    addToRunTimeSelectionTable(shallowFluxScheme, Rusanov, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shallowFluxSchemes::Rusanov::Rusanov
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

Foam::shallowFluxSchemes::Rusanov::~Rusanov()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::shallowFluxSchemes::Rusanov::clear()
{
    shallowFluxScheme::clear();
}


void Foam::shallowFluxSchemes::Rusanov::createSavedFields()
{
    shallowFluxScheme::createSavedFields();
}


void Foam::shallowFluxSchemes::Rusanov::calculateFluxes
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

    scalar lambda(max(mag(UvOwn) + cOwn, mag(UvNei) + cNei));

    scalar phiOwn = UvOwn*magSf;
    scalar phiNei = UvNei*magSf;
    scalar hPhiOwn = hOwn*phiOwn;
    scalar hPhiNei = hNei*phiNei;

    phi = 0.5*(phiOwn + phiNei);
    hPhi = 0.5*(hPhiOwn + hPhiNei - lambda*(hNei - hOwn));
    hUPhi =
        0.5
       *(
            (
                UOwn*hPhiOwn + 0.5*magg_*sqr(hOwn)*normal
              + UNei*hPhiNei + 0.5*magg_*sqr(hNei)*normal
            )*magSf
          - lambda*(hNei*UNei - hOwn*UOwn)
        );
}

// ************************************************************************* //
