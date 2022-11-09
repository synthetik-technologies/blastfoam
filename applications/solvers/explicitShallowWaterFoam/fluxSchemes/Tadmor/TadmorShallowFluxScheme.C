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

#include "TadmorShallowFluxScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace shallowFluxSchemes
{
    defineTypeNameAndDebug(Tadmor, 0);
    addToRunTimeSelectionTable(shallowFluxScheme, Tadmor, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shallowFluxSchemes::Tadmor::Tadmor
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

Foam::shallowFluxSchemes::Tadmor::~Tadmor()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::shallowFluxSchemes::Tadmor::clear()
{
    shallowFluxScheme::clear();
}


void Foam::shallowFluxSchemes::Tadmor::createSavedFields()
{
    shallowFluxScheme::createSavedFields();
}


void Foam::shallowFluxSchemes::Tadmor::calculateFluxes
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

    scalar ap(max(max(UvOwn + cOwn, UvNei + cNei), 0.0));
    scalar am(min(min(UvOwn - cOwn, UvNei - cNei), 0.0));

    scalar amax(max(mag(ap), mag(am)));
    scalar a(-0.5*amax);

    UvOwn *= 0.5;
    UvNei *= 0.5;

    scalar aUvOwn(UvOwn - a);
    scalar aUvNei(UvNei + a);

    phi = magSf*(aUvOwn + aUvNei);

    hPhi = magSf*(aUvOwn*hOwn + aUvNei*hNei);

    hUPhi =
        magSf*(aUvOwn*hOwn*UOwn + aUvNei*hNei*UNei)
      + 0.5*magg_*sqr(0.5*(hOwn + hNei))*Sf;
}


// ************************************************************************* //
