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
    const volScalarField& h0,
    const volVectorField& U,
    const volVectorField& hU
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

    if (!th0Own_.valid())
    {
        autoPtr<ReconstructionScheme<scalar>> h0Limiter
        (
            ReconstructionScheme<scalar>::New(h0, "h0")
        );
        h0Limiter->interpolateOwnNei(th0Own_, th0Nei_);
    }
    const surfaceScalarField& h0Own = th0Own_();
    const surfaceScalarField& h0Nei = th0Nei_();

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
                h0Own[facei], h0Nei[facei],
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
        const vectorField& pSf = mesh_.Sf().boundaryField()[patchi];
        if (bphi[patchi].coupled())
        {
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
                        h0Own.boundaryField()[patchi][facei],
                        h0Nei.boundaryField()[patchi][facei],
                        UOwn.boundaryField()[patchi][facei],
                        UNei.boundaryField()[patchi][facei],
                        pSf[facei],
                        pphi[facei],
                        phPhi[facei],
                        phUPhi[facei],
                        facei, patchi
                    );
                }
            }
        }
        else
        {
            pphi = U.boundaryField()[patchi] & pSf;
            phPhi = hU.boundaryField()[patchi] & pSf;
            phUPhi =
                phPhi*U.boundaryField()[patchi]
              + 0.5*magg_*sqr(h.boundaryField()[patchi])*pSf;
        }
    }
    postUpdate();
}


Foam::tmp<Foam::volVectorField> Foam::shallowFluxScheme::ghGradH0
(
    const dimensionedVector& g,
    const volScalarField& h,
    const volScalarField& h0
) const
{
    return mag(g)*h*fvc::grad(h0);
}
// ************************************************************************* //
