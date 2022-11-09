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

#include "KurganovShallowFluxScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace shallowFluxSchemes
{
    defineTypeNameAndDebug(Kurganov, 0);
    addToRunTimeSelectionTable(shallowFluxScheme, Kurganov, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shallowFluxSchemes::Kurganov::Kurganov
(
    surfaceScalarField& phi,
    surfaceScalarField& hPhi,
    surfaceVectorField& hUPhi,
    const dimensionedVector& g
)
:
    Tadmor(phi, hPhi, hUPhi, g)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::shallowFluxSchemes::Kurganov::~Kurganov()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::shallowFluxSchemes::Kurganov::clear()
{
    Tadmor::clear();
}


void Foam::shallowFluxSchemes::Kurganov::createSavedFields()
{
    Tadmor::createSavedFields();
}


void Foam::shallowFluxSchemes::Kurganov::calculateFluxes
(
    const scalar& hOwn, const scalar& hNei,
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

    scalar ap(max(max(UvOwn + cOwn, UvNei + cNei), 0.0));
    scalar am(min(min(UvOwn - cOwn, UvNei - cNei), 0.0));

    scalar wOwn(ap/(ap - am));
    scalar aSf(am*wOwn);
    scalar wNei(1.0 - wOwn);

    UvOwn *= wOwn;
    UvNei *= wNei;

    scalar aUvOwn(UvOwn - aSf);
    scalar aUvNei(UvNei + aSf);

    // this->save(facei, patchi, aphivOwn, aPhivOwn_);
    // this->save(facei, patchi, aphivNei, aPhivNei_);

    phi = magSf*(aUvOwn + aUvNei);

    hPhi = magSf*(aUvOwn*hOwn + aUvNei*hNei);

    hUPhi =
        magSf*(aUvOwn*hOwn*UOwn + aUvNei*hNei*UNei)
      + 0.5*magg_*sqr(wOwn*hOwn + wNei*hNei)*Sf;
}


// ************************************************************************* //
