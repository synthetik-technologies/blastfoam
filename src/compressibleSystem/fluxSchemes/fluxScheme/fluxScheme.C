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

#include "fluxScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluxScheme, 0);
    defineRunTimeSelectionTable(fluxScheme, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxScheme::fluxScheme(const fvMesh& mesh)
:
    own_
    (
        IOobject
        (
            "own",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        1.0
    ),
    nei_
    (
        IOobject
        (
            "nei",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        -1.0
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxScheme::~fluxScheme()
{}

void Foam::fluxScheme::update
(
    const UPtrList<volScalarField>& alphas,
    const UPtrList<volScalarField>& rhos,
    const volVectorField& U,
    const volScalarField& e,
    const volScalarField& p,
    const volScalarField& c,
    surfaceScalarField& phi,
    UPtrList<surfaceScalarField>& alphaPhis,
    UPtrList<surfaceScalarField>& alphaRhoPhis,
    surfaceScalarField& rhoPhi,
    surfaceVectorField& rhoUPhi,
    surfaceScalarField& rhoEPhi
)
{
    const fvMesh& mesh = U.mesh();

    // Interpolate fields
    PtrList<surfaceScalarField> alphasOwn(alphas.size());
    PtrList<surfaceScalarField> alphasNei(alphas.size());

    PtrList<surfaceScalarField> rhosOwn(alphas.size());
    PtrList<surfaceScalarField> rhosNei(alphas.size());
    surfaceScalarField rhoOwn
    (
        IOobject
        (
            "rhoOwn",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimDensity, 0.0)
    );
    surfaceScalarField rhoNei
    (
        IOobject
        (
            "rhoNei",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimDensity, 0.0)
    );

    forAll(alphas, phasei)
    {
        alphasOwn.set
        (
            phasei,
            fvc::interpolate(alphas[phasei], own_, "reconstruct(alpha)")
        );
        alphasNei.set
        (
            phasei,
            fvc::interpolate(alphas[phasei], nei_, "reconstruct(alpha)")
        );
        rhosOwn.set
        (
            phasei,
            fvc::interpolate(rhos[phasei], own_, "reconstruct(rho)")
        );
        rhosNei.set
        (
            phasei,
            fvc::interpolate(rhos[phasei], nei_, "reconstruct(rho)")
        );
        rhoOwn += alphasOwn[phasei]*rhosOwn[phasei];
        rhoNei += alphasNei[phasei]*rhosNei[phasei];
    }

    surfaceVectorField UOwn(fvc::interpolate(U, own_, "reconstruct(U)"));
    surfaceVectorField UNei(fvc::interpolate(U, nei_, "reconstruct(U)"));

    surfaceScalarField pOwn(fvc::interpolate(p, own_, "reconstruct(p)"));
    surfaceScalarField pNei(fvc::interpolate(p, nei_, "reconstruct(p)"));

    surfaceScalarField cOwn(fvc::interpolate(c, own_, "reconstruct(c)"));
    surfaceScalarField cNei(fvc::interpolate(c, nei_, "reconstruct(c)"));

    surfaceScalarField eOwn(fvc::interpolate(e, own_, "reconstruct(e)"));
    surfaceScalarField eNei(fvc::interpolate(e, nei_, "reconstruct(e)"));

    forAll(UOwn, facei)
    {
        scalarList alphasiOwn(alphas.size());
        scalarList alphasiNei(alphas.size());
        scalarList rhosiOwn(alphas.size());
        scalarList rhosiNei(alphas.size());

        scalarList alphaPhisi(alphas.size());
        scalarList alphaRhoPhisi(alphas.size());

        forAll(alphas, phasei)
        {
            alphasiOwn[phasei] = alphasOwn[phasei][facei];
            alphasiNei[phasei] = alphasNei[phasei][facei];
            rhosiOwn[phasei] = rhosOwn[phasei][facei];
            rhosiNei[phasei] = rhosNei[phasei][facei];
        }
        calculateFluxes
        (
            alphasiOwn, alphasiNei,
            rhosiOwn, rhosiNei,
            rhoOwn[facei], rhoNei[facei],
            UOwn[facei], UNei[facei],
            eOwn[facei], eNei[facei],
            pOwn[facei], pNei[facei],
            cOwn[facei], cNei[facei],
            mesh.Sf()[facei],
            phi[facei],
            alphaPhisi,
            alphaRhoPhisi,
            rhoUPhi[facei],
            rhoEPhi[facei]
        );

        rhoPhi[facei] = 0.0;
        forAll(alphas, phasei)
        {
            alphaPhis[phasei][facei] = alphaPhisi[phasei];
            alphaRhoPhis[phasei][facei] = alphaRhoPhisi[phasei];
            rhoPhi[facei] += alphaRhoPhisi[phasei];
        }
    }

    forAll(U.boundaryField(), patchi)
    {
        forAll(U.boundaryField()[patchi], facei)
        {
            scalarList alphasiOwn(alphas.size());
            scalarList alphasiNei(alphas.size());
            scalarList rhosiOwn(alphas.size());
            scalarList rhosiNei(alphas.size());

            scalarList alphaPhisi(alphas.size());
            scalarList alphaRhoPhisi(alphas.size());

            forAll(alphas, phasei)
            {
                alphasiOwn[phasei] =
                    alphasOwn[phasei].boundaryField()[patchi][facei];
                alphasiNei[phasei] =
                    alphasNei[phasei].boundaryField()[patchi][facei];
                rhosiOwn[phasei] =
                    rhosOwn[phasei].boundaryField()[patchi][facei];
                rhosiNei[phasei] =
                    rhosNei[phasei].boundaryField()[patchi][facei];
            }
            calculateFluxes
            (
                alphasiOwn, alphasiNei,
                rhosiOwn, rhosiNei,
                rhoOwn.boundaryField()[patchi][facei],
                rhoNei.boundaryField()[patchi][facei],
                UOwn.boundaryField()[patchi][facei],
                UNei.boundaryField()[patchi][facei],
                eOwn.boundaryField()[patchi][facei],
                eNei.boundaryField()[patchi][facei],
                pOwn.boundaryField()[patchi][facei],
                pNei.boundaryField()[patchi][facei],
                cOwn.boundaryField()[patchi][facei],
                cNei.boundaryField()[patchi][facei],
                mesh.Sf().boundaryField()[patchi][facei],
                phi.boundaryFieldRef()[patchi][facei],
                alphaPhisi,
                alphaRhoPhisi,
                rhoUPhi.boundaryFieldRef()[patchi][facei],
                rhoEPhi.boundaryFieldRef()[patchi][facei]
            );

            rhoPhi.boundaryFieldRef()[patchi][facei] = 0.0;
            forAll(alphas, phasei)
            {
                alphaPhis[phasei].boundaryFieldRef()[patchi][facei] =
                    alphaPhisi[phasei];
                alphaRhoPhis[phasei].boundaryFieldRef()[patchi][facei] =
                    alphaRhoPhisi[phasei];
                rhoPhi.boundaryFieldRef()[patchi][facei] +=
                    alphaRhoPhisi[phasei];
            }
        }
    }
}


void Foam::fluxScheme::update
(
    const volScalarField& alpha,
    const volScalarField& rho1,
    const volScalarField& rho2,
    const volVectorField& U,
    const volScalarField& e,
    const volScalarField& p,
    const volScalarField& c,
    surfaceScalarField& phi,
    surfaceScalarField& alphaPhi,
    surfaceScalarField& alphaRhoPhi1,
    surfaceScalarField& alphaRhoPhi2,
    surfaceScalarField& rhoPhi,
    surfaceVectorField& rhoUPhi,
    surfaceScalarField& rhoEPhi
)
{
    const fvMesh& mesh = U.mesh();

    // Interpolate fields
    surfaceScalarField alphaOwn
    (
        fvc::interpolate(alpha, own_, "reconstruct(alpha)")
    );
    surfaceScalarField alphaNei
    (
        fvc::interpolate(alpha, nei_, "reconstruct(alpha)")
    );

    surfaceScalarField rho1Own
    (
        fvc::interpolate(rho1, own_, "reconstruct(rho)")
    );

    surfaceScalarField rho1Nei
    (
        fvc::interpolate(rho1, nei_, "reconstruct(rho)")
    );

    surfaceScalarField rho2Own
    (
        fvc::interpolate(rho2, own_, "reconstruct(rho)")
    );

    surfaceScalarField rho2Nei
    (
        fvc::interpolate(rho2, nei_, "reconstruct(rho)")
    );
    surfaceScalarField rhoOwn(alphaOwn*rho1Own + (1.0 - alphaOwn)*rho2Own);
    surfaceScalarField rhoNei(alphaNei*rho1Nei + (1.0 - alphaNei)*rho2Nei);

    surfaceVectorField UOwn(fvc::interpolate(U, own_, "reconstruct(U)"));
    surfaceVectorField UNei(fvc::interpolate(U, nei_, "reconstruct(U)"));

    surfaceScalarField pOwn(fvc::interpolate(p, own_, "reconstruct(p)"));
    surfaceScalarField pNei(fvc::interpolate(p, nei_, "reconstruct(p)"));

    surfaceScalarField cOwn(fvc::interpolate(c, own_, "reconstruct(c)"));
    surfaceScalarField cNei(fvc::interpolate(c, nei_, "reconstruct(c)"));

    surfaceScalarField eOwn(fvc::interpolate(e, own_, "reconstruct(e)"));
    surfaceScalarField eNei(fvc::interpolate(e, nei_, "reconstruct(e)"));

    forAll(UOwn, facei)
    {
        scalarList alphaPhisi(2);
        scalarList alphaRhoPhisi(2);
        calculateFluxes
        (
            {alphaOwn[facei], 1.0 - alphaOwn[facei]},
            {alphaNei[facei], 1.0 - alphaNei[facei]},
            {rho1Own[facei], rho2Own[facei]},
            {rho1Nei[facei], rho2Nei[facei]},
            rhoOwn[facei], rhoNei[facei],
            UOwn[facei], UNei[facei],
            eOwn[facei], eNei[facei],
            pOwn[facei], pNei[facei],
            cOwn[facei], cNei[facei],
            mesh.Sf()[facei],
            phi[facei],
            alphaPhisi,
            alphaRhoPhisi,
            rhoUPhi[facei],
            rhoEPhi[facei]
        );

        alphaPhi[facei] = alphaPhisi[0];
        alphaRhoPhi1[facei] = alphaRhoPhisi[0];
        alphaRhoPhi2[facei] = alphaRhoPhisi[1];

        rhoPhi[facei] = alphaRhoPhi1[facei] + alphaRhoPhi2[facei];
    }

    forAll(U.boundaryField(), patchi)
    {
        forAll(U.boundaryField()[patchi], facei)
        {
            scalarList alphaPhisi(2);
            scalarList alphaRhoPhisi(2);

            calculateFluxes
            (
                {
                    alphaOwn.boundaryField()[patchi][facei],
                    1.0 - alphaOwn.boundaryField()[patchi][facei]
                },
                {
                    alphaNei.boundaryField()[patchi][facei],
                    1.0 - alphaNei.boundaryField()[patchi][facei]
                },
                {
                    rho1Own.boundaryField()[patchi][facei],
                    rho2Own.boundaryField()[patchi][facei]
                },
                {
                    rho1Nei.boundaryField()[patchi][facei],
                    rho2Nei.boundaryField()[patchi][facei]
                },
                rhoOwn.boundaryField()[patchi][facei],
                rhoNei.boundaryField()[patchi][facei],
                UOwn.boundaryField()[patchi][facei],
                UNei.boundaryField()[patchi][facei],
                eOwn.boundaryField()[patchi][facei],
                eNei.boundaryField()[patchi][facei],
                pOwn.boundaryField()[patchi][facei],
                pNei.boundaryField()[patchi][facei],
                cOwn.boundaryField()[patchi][facei],
                cNei.boundaryField()[patchi][facei],
                mesh.Sf().boundaryField()[patchi][facei],
                phi.boundaryFieldRef()[patchi][facei],
                alphaPhisi,
                alphaRhoPhisi,
                rhoUPhi.boundaryFieldRef()[patchi][facei],
                rhoEPhi.boundaryFieldRef()[patchi][facei]
            );

            alphaPhi.boundaryFieldRef()[patchi][facei] = alphaPhisi[0];
            alphaRhoPhi1.boundaryFieldRef()[patchi][facei] = alphaRhoPhisi[0];
            alphaRhoPhi2.boundaryFieldRef()[patchi][facei] = alphaRhoPhisi[1];

            rhoPhi.boundaryFieldRef()[patchi][facei] =
                alphaRhoPhi1.boundaryField()[patchi][facei]
              + alphaRhoPhi2.boundaryField()[patchi][facei];
        }
    }
}

// ************************************************************************* //
