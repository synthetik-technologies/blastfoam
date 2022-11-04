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

#include "shallowFluxScheme.H"
#include "ReconstructionScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(shallowFluxScheme, 0);
    defineRunTimeSelectionTable(shallowFluxScheme, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shallowFluxScheme::shallowFluxScheme
(
    surfaceScalarField& phi,
    surfaceScalarField& hPhi,
    surfaceVectorField& hUPhi,
    const dimensionedVector& g
)
:
    fluxSchemeBase(phi),
    dict_(mesh_.schemesDict()),
    phi_(phi),
    hPhi_(hPhi),
    hUPhi_(hUPhi),
    g_(g.value()),
    magg_(mag(g_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::shallowFluxScheme::~shallowFluxScheme()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::shallowFluxScheme::clear()
{}

void Foam::shallowFluxScheme::createSavedFields()
{}


void Foam::shallowFluxScheme::update
(
    const volScalarField& h,
    const volVectorField& U
)
{
    createSavedFields();

    autoPtr<ReconstructionScheme<scalar>> hLimiter
    (
        ReconstructionScheme<scalar>::New(h, "h")
    );
    autoPtr<ReconstructionScheme<vector>> ULimiter
    (
        ReconstructionScheme<vector>::New(U, "U")
    );

    tmp<surfaceScalarField> thOwn, thNei;
    hLimiter->interpolateOwnNei(thOwn, thNei);
    const surfaceScalarField& hOwn = thOwn();
    const surfaceScalarField& hNei = thNei();

    tmp<surfaceVectorField> tUOwn, tUNei;
    ULimiter->interpolateOwnNei(tUOwn, tUNei);
    const surfaceVectorField& UOwn = tUOwn();
    const surfaceVectorField& UNei = tUNei();

    scalarField& phiI = phi_.primitiveFieldRef();
    scalarField& hPhiI = hPhi_.primitiveFieldRef();
    vectorField& hUPhiI = hUPhi_.primitiveFieldRef();
    // preUpdate(p);
    forAll(UOwn, facei)
    {
        if (max(hOwn[facei], hNei[facei]) > small)
        {
            calculateFluxes
            (
                hOwn[facei], hNei[facei],
                UOwn[facei], UNei[facei],
                mesh_.Sf()[facei],
                phiI[facei],
                hPhiI[facei],
                hUPhiI[facei],
                facei
            );
        }
    }

    surfaceScalarField::Boundary& bphi = phi_.boundaryFieldRef();
    surfaceScalarField::Boundary& bhPhi = hPhi_.boundaryFieldRef();
    surfaceVectorField::Boundary& bhUPhi = hUPhi_.boundaryFieldRef();
    forAll(bphi, patchi)
    {
        scalarField& pphi = bphi[patchi];
        scalarField& phPhi = bhPhi[patchi];
        vectorField& phUPhi = bhUPhi[patchi];
        forAll(U.boundaryField()[patchi], facei)
        {
            if
            (
                max
                (
                    hOwn.boundaryField()[patchi][facei],
                    hNei.boundaryField()[patchi][facei]
                ) > small
            )
            {
                calculateFluxes
                (
                    hOwn.boundaryField()[patchi][facei],
                    hNei.boundaryField()[patchi][facei],
                    UOwn.boundaryField()[patchi][facei],
                    UNei.boundaryField()[patchi][facei],
                    mesh_.Sf().boundaryField()[patchi][facei],
                    pphi[facei],
                    phPhi[facei],
                    phUPhi[facei],
                    facei, patchi
                );
            }
        }
    }
    postUpdate();
}

// ************************************************************************* //
