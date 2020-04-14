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
    regIOobject
    (
        IOobject
        (
            "fluxScheme",
            mesh.time().timeName(),
            mesh
        )
    ),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxScheme::~fluxScheme()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluxScheme::clear()
{
    own_.clear();
    nei_.clear();
    Uf_.clear();
    rhoOwn_.clear();
    rhoNei_.clear();
}

void Foam::fluxScheme::createSavedFields()
{
    if (own_.valid())
    {
        return;
    }
    own_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "fluxScheme::own",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("1", dimless, 1.0)
        )
    );
    nei_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "fluxScheme::nei",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("-1", dimless, -1.0)
        )
    );
    Uf_ = tmp<surfaceVectorField>
    (
        new surfaceVectorField
        (
            IOobject
            (
                "fluxScheme::Uf",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedVector("0", dimVelocity, Zero)
        )
    );
}

Foam::tmp<Foam::surfaceVectorField> Foam::fluxScheme::Uf() const
{
    if (Uf_.valid())
    {
        return Uf_;
    }
    return tmp<surfaceVectorField>
    (
        new surfaceVectorField
        (
            IOobject
            (
                "fluxScheme::Uf",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedVector("0", dimVelocity, Zero)
        )
    );
}


void Foam::fluxScheme::update
(
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& e,
    const volScalarField& p,
    const volScalarField& c,
    surfaceScalarField& phi,
    surfaceScalarField& rhoPhi,
    surfaceVectorField& rhoUPhi,
    surfaceScalarField& rhoEPhi
)
{
    createSavedFields();

    rhoOwn_ = fvc::interpolate(rho, own_(), scheme("rho"));
    rhoNei_ = fvc::interpolate(rho, nei_(), scheme("rho"));

    surfaceVectorField UOwn(fvc::interpolate(U, own_(), scheme("U")));
    surfaceVectorField UNei(fvc::interpolate(U, nei_(), scheme("U")));

    surfaceScalarField pOwn(fvc::interpolate(p, own_(), scheme("p")));
    surfaceScalarField pNei(fvc::interpolate(p, nei_(), scheme("p")));

    surfaceScalarField cOwn(fvc::interpolate(c, own_(), scheme("c")));
    surfaceScalarField cNei(fvc::interpolate(c, nei_(), scheme("c")));

    surfaceScalarField eOwn(fvc::interpolate(e, own_(), scheme("e")));
    surfaceScalarField eNei(fvc::interpolate(e, nei_(), scheme("e")));

    preUpdate(p);
    forAll(UOwn, facei)
    {

        calculateFluxes
        (
            rhoOwn_()[facei], rhoNei_()[facei],
            UOwn[facei], UNei[facei],
            eOwn[facei], eNei[facei],
            pOwn[facei], pNei[facei],
            cOwn[facei], cNei[facei],
            mesh_.Sf()[facei],
            phi[facei],
            rhoPhi[facei],
            rhoUPhi[facei],
            rhoEPhi[facei],
            facei
        );
    }

    forAll(U.boundaryField(), patchi)
    {
        forAll(U.boundaryField()[patchi], facei)
        {

            calculateFluxes
            (
                rhoOwn_().boundaryField()[patchi][facei],
                rhoNei_().boundaryField()[patchi][facei],
                UOwn.boundaryField()[patchi][facei],
                UNei.boundaryField()[patchi][facei],
                eOwn.boundaryField()[patchi][facei],
                eNei.boundaryField()[patchi][facei],
                pOwn.boundaryField()[patchi][facei],
                pNei.boundaryField()[patchi][facei],
                cOwn.boundaryField()[patchi][facei],
                cNei.boundaryField()[patchi][facei],
                mesh_.Sf().boundaryField()[patchi][facei],
                phi.boundaryFieldRef()[patchi][facei],
                rhoPhi.boundaryFieldRef()[patchi][facei],
                rhoUPhi.boundaryFieldRef()[patchi][facei],
                rhoEPhi.boundaryFieldRef()[patchi][facei],
                facei, patchi
            );
        }
    }
    postUpdate();
}

void Foam::fluxScheme::update
(
    const PtrList<volScalarField>& alphas,
    const UPtrList<volScalarField>& rhos,
    const volVectorField& U,
    const volScalarField& e,
    const volScalarField& p,
    const volScalarField& c,
    surfaceScalarField& phi,
    PtrList<surfaceScalarField>& alphaPhis,
    PtrList<surfaceScalarField>& alphaRhoPhis,
    surfaceScalarField& rhoPhi,
    surfaceVectorField& rhoUPhi,
    surfaceScalarField& rhoEPhi
)
{
    createSavedFields();

    // Interpolate fields
    PtrList<surfaceScalarField> alphasOwn(alphas.size());
    PtrList<surfaceScalarField> alphasNei(alphas.size());

    PtrList<surfaceScalarField> rhosOwn(alphas.size());
    PtrList<surfaceScalarField> rhosNei(alphas.size());
    rhoOwn_ =
    (
        new surfaceScalarField
        (
            IOobject
            (
                "rhoOwn",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimDensity, 0.0)
        )
    );
    rhoNei_ =
    (
        new surfaceScalarField
        (
            IOobject
            (
                "rhoNei",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimDensity, 0.0)
        )
    );

    forAll(alphas, phasei)
    {
        alphasOwn.set
        (
            phasei,
            fvc::interpolate(alphas[phasei], own_(), scheme("alpha"))
        );
        alphasNei.set
        (
            phasei,
            fvc::interpolate(alphas[phasei], nei_(), scheme("alpha"))
        );
        rhosOwn.set
        (
            phasei,
            fvc::interpolate(rhos[phasei], own_(), scheme("rho"))
        );
        rhosNei.set
        (
            phasei,
            fvc::interpolate(rhos[phasei], nei_(), scheme("rho"))
        );
        rhoOwn_.ref() += alphasOwn[phasei]*rhosOwn[phasei];
        rhoNei_.ref() += alphasNei[phasei]*rhosNei[phasei];
    }

    surfaceVectorField UOwn(fvc::interpolate(U, own_(), scheme("U")));
    surfaceVectorField UNei(fvc::interpolate(U, nei_(), scheme("U")));

    surfaceScalarField pOwn(fvc::interpolate(p, own_(), scheme("p")));
    surfaceScalarField pNei(fvc::interpolate(p, nei_(), scheme("p")));

    surfaceScalarField cOwn(fvc::interpolate(c, own_(), scheme("c")));
    surfaceScalarField cNei(fvc::interpolate(c, nei_(), scheme("c")));

    surfaceScalarField eOwn(fvc::interpolate(e, own_(), scheme("e")));
    surfaceScalarField eNei(fvc::interpolate(e, nei_(), scheme("e")));

    preUpdate(p);
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
            rhoOwn_()[facei], rhoNei_()[facei],
            UOwn[facei], UNei[facei],
            eOwn[facei], eNei[facei],
            pOwn[facei], pNei[facei],
            cOwn[facei], cNei[facei],
            mesh_.Sf()[facei],
            phi[facei],
            alphaPhisi,
            alphaRhoPhisi,
            rhoUPhi[facei],
            rhoEPhi[facei],
            facei
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
                rhoOwn_().boundaryField()[patchi][facei],
                rhoNei_().boundaryField()[patchi][facei],
                UOwn.boundaryField()[patchi][facei],
                UNei.boundaryField()[patchi][facei],
                eOwn.boundaryField()[patchi][facei],
                eNei.boundaryField()[patchi][facei],
                pOwn.boundaryField()[patchi][facei],
                pNei.boundaryField()[patchi][facei],
                cOwn.boundaryField()[patchi][facei],
                cNei.boundaryField()[patchi][facei],
                mesh_.Sf().boundaryField()[patchi][facei],
                phi.boundaryFieldRef()[patchi][facei],
                alphaPhisi,
                alphaRhoPhisi,
                rhoUPhi.boundaryFieldRef()[patchi][facei],
                rhoEPhi.boundaryFieldRef()[patchi][facei],
                facei, patchi
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
    postUpdate();
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
    createSavedFields();

    // Interpolate fields
    surfaceScalarField alphaOwn
    (
        fvc::interpolate(alpha, own_(), scheme("alpha"))
    );
    surfaceScalarField alphaNei
    (
        fvc::interpolate(alpha, nei_(), scheme("alpha"))
    );

    surfaceScalarField rho1Own
    (
        fvc::interpolate(rho1, own_(), scheme("rho"))
    );

    surfaceScalarField rho1Nei
    (
        fvc::interpolate(rho1, nei_(), scheme("rho"))
    );

    surfaceScalarField rho2Own
    (
        fvc::interpolate(rho2, own_(), scheme("rho"))
    );

    surfaceScalarField rho2Nei
    (
        fvc::interpolate(rho2, nei_(), scheme("rho"))
    );
    rhoOwn_ = (alphaOwn*rho1Own + (1.0 - alphaOwn)*rho2Own);
    rhoNei_ = (alphaNei*rho1Nei + (1.0 - alphaNei)*rho2Nei);

    surfaceVectorField UOwn(fvc::interpolate(U, own_(), scheme("U")));
    surfaceVectorField UNei(fvc::interpolate(U, nei_(), scheme("U")));

    surfaceScalarField pOwn(fvc::interpolate(p, own_(), scheme("p")));
    surfaceScalarField pNei(fvc::interpolate(p, nei_(), scheme("p")));

    surfaceScalarField cOwn(fvc::interpolate(c, own_(), scheme("c")));
    surfaceScalarField cNei(fvc::interpolate(c, nei_(), scheme("c")));

    surfaceScalarField eOwn(fvc::interpolate(e, own_(), scheme("e")));
    surfaceScalarField eNei(fvc::interpolate(e, nei_(), scheme("e")));

    preUpdate(p);
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
            rhoOwn_()[facei], rhoNei_()[facei],
            UOwn[facei], UNei[facei],
            eOwn[facei], eNei[facei],
            pOwn[facei], pNei[facei],
            cOwn[facei], cNei[facei],
            mesh_.Sf()[facei],
            phi[facei],
            alphaPhisi,
            alphaRhoPhisi,
            rhoUPhi[facei],
            rhoEPhi[facei],
            facei
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
                rhoOwn_().boundaryField()[patchi][facei],
                rhoNei_().boundaryField()[patchi][facei],
                UOwn.boundaryField()[patchi][facei],
                UNei.boundaryField()[patchi][facei],
                eOwn.boundaryField()[patchi][facei],
                eNei.boundaryField()[patchi][facei],
                pOwn.boundaryField()[patchi][facei],
                pNei.boundaryField()[patchi][facei],
                cOwn.boundaryField()[patchi][facei],
                cNei.boundaryField()[patchi][facei],
                mesh_.Sf().boundaryField()[patchi][facei],
                phi.boundaryFieldRef()[patchi][facei],
                alphaPhisi,
                alphaRhoPhisi,
                rhoUPhi.boundaryFieldRef()[patchi][facei],
                rhoEPhi.boundaryFieldRef()[patchi][facei],
                facei, patchi
            );

            alphaPhi.boundaryFieldRef()[patchi][facei] = alphaPhisi[0];
            alphaRhoPhi1.boundaryFieldRef()[patchi][facei] = alphaRhoPhisi[0];
            alphaRhoPhi2.boundaryFieldRef()[patchi][facei] = alphaRhoPhisi[1];

            rhoPhi.boundaryFieldRef()[patchi][facei] =
                alphaRhoPhi1.boundaryField()[patchi][facei]
              + alphaRhoPhi2.boundaryField()[patchi][facei];
        }
    }
    postUpdate();
}


Foam::tmp<Foam::surfaceScalarField> Foam::fluxScheme::energyFlux
(
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& e,
    const volScalarField& p
) const
{
    tmp<surfaceScalarField> rhoOwn;
    tmp<surfaceScalarField> rhoNei;
    if (rho.name() == "rho")
    {
        rhoOwn = tmp<surfaceScalarField>(new surfaceScalarField(rhoOwn_()));
        rhoNei = tmp<surfaceScalarField>(new surfaceScalarField(rhoNei_()));
    }
    else
    {
        rhoOwn = fvc::interpolate(rho, own_(), scheme("rho"));
        rhoNei = fvc::interpolate(rho, nei_(), scheme("rho"));
    }

    surfaceVectorField UOwn
    (
        fvc::interpolate(U, own_(), scheme("U"))
    );
    surfaceVectorField UNei
    (
        fvc::interpolate(U, nei_(), scheme("U"))
    );

    surfaceScalarField eOwn
    (
        fvc::interpolate(e, own_(), scheme("e"))
    );
    surfaceScalarField eNei
    (
        fvc::interpolate(e, nei_(), scheme("e"))
    );

    surfaceScalarField pOwn
    (
        fvc::interpolate(p, own_(), scheme("p"))
    );
    surfaceScalarField pNei
    (
        fvc::interpolate(p, nei_(), scheme("p"))
    );

    tmp<surfaceScalarField> tmpPhi
    (
        new surfaceScalarField
        (
            IOobject
            (
                e.name() + "Phi",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar
            (
                "0",
                rho.dimensions()*e.dimensions()*dimVelocity*dimArea,
                0.0
            )
        )
    );
    surfaceScalarField& phi = tmpPhi.ref();

    forAll(eOwn, facei)
    {
        phi[facei] = energyFlux
        (
            rhoOwn()[facei], rhoNei()[facei],
            UOwn[facei], UNei[facei],
            eOwn[facei], eNei[facei],
            pOwn[facei], pNei[facei],
            facei
        );
    }

    forAll(e.boundaryField(), patchi)
    {
        forAll(e.boundaryField()[patchi], facei)
        {
            phi.boundaryFieldRef()[patchi][facei] =
                energyFlux
                (
                    rhoOwn().boundaryField()[patchi][facei],
                    rhoNei().boundaryField()[patchi][facei],
                    UOwn.boundaryField()[patchi][facei],
                    UNei.boundaryField()[patchi][facei],
                    eOwn.boundaryField()[patchi][facei],
                    eNei.boundaryField()[patchi][facei],
                    pOwn.boundaryField()[patchi][facei],
                    pNei.boundaryField()[patchi][facei],
                    facei, patchi
                );
        }
    }
    return tmpPhi;
}


bool Foam::fluxScheme::writeData(Ostream& os) const
{
    return os.good();
}

// ************************************************************************* //
