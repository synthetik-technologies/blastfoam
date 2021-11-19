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
#include "MUSCLReconstructionScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluxScheme, 0);
    defineRunTimeSelectionTable(fluxScheme, singlePhase);
    defineRunTimeSelectionTable(fluxScheme, multiphase);
    defineRunTimeSelectionTable(fluxScheme, interface);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxScheme::fluxScheme(const fvMesh& mesh)
:
    fluxSchemeBase(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxScheme::~fluxScheme()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluxScheme::clear()
{
    Uf_.clear();
}

void Foam::fluxScheme::createSavedFields()
{
    if (Uf_.valid())
    {
        return;
    }
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
    return surfaceVectorField::New
    (
        "fluxScheme::Uf",
        mesh_,
        dimensionedVector("0", dimVelocity, Zero)
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::fluxScheme::AD(const volScalarField& alpha) const
{
    return mag((this->Uf_() & alpha.mesh().Sf())/alpha.mesh().magSf());
}


Foam::tmp<Foam::surfaceScalarField>
Foam::fluxScheme::snGradAlpha
(
    const volScalarField& alpha
) const
{
    const fvMesh& mesh = alpha.mesh();
    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();
    surfaceVectorField n(Sf/magSf);

    surfaceScalarField alphaf
    (
        surfaceInterpolationScheme<scalar>::New
        (
            mesh,
            IStringStream("linear vanLeer")()
        )->interpolate(alpha)
    );
    volVectorField gradAlpha(fvc::grad(alphaf));
    tmp<surfaceScalarField> tsnGradAlpha
    (
        surfaceScalarField::New
        (
            IOobject::groupName("snGradAlpha", alpha.group()),
            fvc::snGrad(alpha)
        )
    );
    surfaceScalarField& snGradAlpha = tsnGradAlpha.ref();
    forAll(alphaf, facei)
    {
        if
        (
            ((gradAlpha[own[facei]] & n[facei])*snGradAlpha[facei]) > 0
         && mag(gradAlpha[own[facei]] & n[facei])
          < mag(snGradAlpha[facei])
        )
        {
            alphaf[facei] = alpha[nei[facei]];
        }
        else if
        (
            mag
            (
                mag(gradAlpha[own[facei]] & n[facei])
              - mag(snGradAlpha[facei])
            ) < small
        )
        {
            alphaf[facei] = 0.5*(alpha[own[facei]] + alpha[nei[facei]]);
        }
        else
        {
            alphaf[facei] = alpha[own[facei]];
        }
    }
    gradAlpha = fvc::grad(alphaf);

    forAll(snGradAlpha, facei)
    {
        snGradAlpha[facei] =
            minMagSqrOp<vector>()
            (
                gradAlpha[own[facei]],
                gradAlpha[nei[facei]]
            ) & Sf[facei];
    }
    surfaceScalarField::Boundary& bsnGradAlpha
    (
        snGradAlpha.boundaryFieldRef()
    );
    forAll(bsnGradAlpha, patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        const fvPatchField<vector>& pgradAlpha
        (
            gradAlpha.boundaryField()[patchi]
        );
        const vectorField& pSf(Sf.boundaryField()[patchi]);
        scalarField& psnGradAlpha(bsnGradAlpha[patchi]);

        if (patch.coupled())
        {
            vectorField gradAlphaOwn(pgradAlpha.patchInternalField());
            vectorField gradAlphaNei(pgradAlpha.patchNeighbourField());
            forAll(gradAlphaOwn, facei)
            {
                snGradAlpha[facei] =
                    minMagSqrOp<vector>()
                    (
                        gradAlphaOwn[facei],
                        gradAlphaNei[facei]
                    ) & pSf[facei];
            }
        }
        else
        {
            psnGradAlpha = pgradAlpha & pSf;
        }
    }
    snGradAlpha.dimensions().reset(snGradAlpha.dimensions()*dimVolume);
    return tsnGradAlpha;
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

    autoPtr<MUSCLReconstructionScheme<scalar>> rhoLimiter
    (
        MUSCLReconstructionScheme<scalar>::New(rho, "rho")
    );
    autoPtr<MUSCLReconstructionScheme<vector>> ULimiter
    (
        MUSCLReconstructionScheme<vector>::New(U, "U")
    );
    autoPtr<MUSCLReconstructionScheme<scalar>> eLimiter
    (
        MUSCLReconstructionScheme<scalar>::New(e, "e")
    );
    autoPtr<MUSCLReconstructionScheme<scalar>> pLimiter
    (
        MUSCLReconstructionScheme<scalar>::New(p, "p")
    );
    autoPtr<MUSCLReconstructionScheme<scalar>> cLimiter
    (
        MUSCLReconstructionScheme<scalar>::New(c, "speedOfSound")
    );

    tmp<surfaceScalarField> trhoOwn(rhoLimiter->interpolateOwn());
    tmp<surfaceScalarField> trhoNei(rhoLimiter->interpolateNei());
    const surfaceScalarField& rhoOwn = trhoOwn();
    const surfaceScalarField& rhoNei = trhoNei();

    tmp<surfaceVectorField> tUOwn(ULimiter->interpolateOwn());
    tmp<surfaceVectorField> tUNei(ULimiter->interpolateNei());
    const surfaceVectorField& UOwn = tUOwn();
    const surfaceVectorField& UNei = tUNei();

    tmp<surfaceScalarField> teOwn(eLimiter->interpolateOwn());
    tmp<surfaceScalarField> teNei(eLimiter->interpolateNei());
    const surfaceScalarField& eOwn = teOwn();
    const surfaceScalarField& eNei = teNei();

    tmp<surfaceScalarField> tpOwn(pLimiter->interpolateOwn());
    tmp<surfaceScalarField> tpNei(pLimiter->interpolateNei());
    const surfaceScalarField& pOwn = tpOwn();
    const surfaceScalarField& pNei = tpNei();

    tmp<surfaceScalarField> tcOwn(cLimiter->interpolateOwn());
    tmp<surfaceScalarField> tcNei(cLimiter->interpolateNei());
    const surfaceScalarField& cOwn = tcOwn();
    const surfaceScalarField& cNei = tcNei();


    preUpdate(p);
    forAll(UOwn, facei)
    {

        calculateFluxes
        (
            rhoOwn[facei], rhoNei[facei],
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
        scalarField& pphi = phi.boundaryFieldRef()[patchi];
        scalarField& prhoPhi = rhoPhi.boundaryFieldRef()[patchi];
        vectorField& prhoUPhi = rhoUPhi.boundaryFieldRef()[patchi];
        scalarField& prhoEPhi = rhoEPhi.boundaryFieldRef()[patchi];
        forAll(U.boundaryField()[patchi], facei)
        {
            calculateFluxes
            (
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
                mesh_.Sf().boundaryField()[patchi][facei],
                pphi[facei],
                prhoPhi[facei],
                prhoUPhi[facei],
                prhoEPhi[facei],
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
    tmp<surfaceScalarField> trhoOwn
    (
        surfaceScalarField::New
        (
            "rhoOwn",
            mesh_,
            dimensionedScalar("0", dimDensity, 0.0)
        )
    );
    surfaceScalarField& rhoOwn = trhoOwn.ref();
    tmp<surfaceScalarField> trhoNei
    (
        surfaceScalarField::New
        (
            "rhoNei",
            mesh_,
            dimensionedScalar("0", dimDensity, 0.0)
        )
    );
    surfaceScalarField& rhoNei = trhoNei.ref();

    forAll(alphas, phasei)
    {
        const word phaseName(alphas[phasei].group());
        autoPtr<MUSCLReconstructionScheme<scalar>> alphaLimiter
        (
            MUSCLReconstructionScheme<scalar>::New
            (
                alphas[phasei],
                "alpha",
                phaseName
            )
        );
        autoPtr<MUSCLReconstructionScheme<scalar>> rhoLimiter
        (
            MUSCLReconstructionScheme<scalar>::New
            (
                rhos[phasei],
                "rho",
                phaseName
            )
        );
        alphasOwn.set
        (
            phasei,
            alphaLimiter->interpolateOwn()
        );
        alphasNei.set
        (
            phasei,
            alphaLimiter->interpolateNei()
        );
        rhosOwn.set
        (
            phasei,
            rhoLimiter->interpolateOwn()
        );
        rhosNei.set
        (
            phasei,
            rhoLimiter->interpolateNei()
        );
        rhoOwn += alphasOwn[phasei]*rhosOwn[phasei];
        rhoNei += alphasNei[phasei]*rhosNei[phasei];
    }

    autoPtr<MUSCLReconstructionScheme<vector>> ULimiter
    (
        MUSCLReconstructionScheme<vector>::New(U, "U")
    );
    autoPtr<MUSCLReconstructionScheme<scalar>> eLimiter
    (
        MUSCLReconstructionScheme<scalar>::New(e, "e")
    );
    autoPtr<MUSCLReconstructionScheme<scalar>> pLimiter
    (
        MUSCLReconstructionScheme<scalar>::New(p, "p")
    );
    autoPtr<MUSCLReconstructionScheme<scalar>> cLimiter
    (
        MUSCLReconstructionScheme<scalar>::New(c, "speedOfSound")
    );

    tmp<surfaceVectorField> tUOwn(ULimiter->interpolateOwn());
    tmp<surfaceVectorField> tUNei(ULimiter->interpolateNei());
    const surfaceVectorField& UOwn = tUOwn();
    const surfaceVectorField& UNei = tUNei();

    tmp<surfaceScalarField> teOwn(eLimiter->interpolateOwn());
    tmp<surfaceScalarField> teNei(eLimiter->interpolateNei());
    const surfaceScalarField& eOwn = teOwn();
    const surfaceScalarField& eNei = teNei();

    tmp<surfaceScalarField> tpOwn(pLimiter->interpolateOwn());
    tmp<surfaceScalarField> tpNei(pLimiter->interpolateNei());
    const surfaceScalarField& pOwn = tpOwn();
    const surfaceScalarField& pNei = tpNei();

    tmp<surfaceScalarField> tcOwn(cLimiter->interpolateOwn());
    tmp<surfaceScalarField> tcNei(cLimiter->interpolateNei());
    const surfaceScalarField& cOwn = tcOwn();
    const surfaceScalarField& cNei = tcNei();

    preUpdate(p);

    // Allocate list for cell and face operations
    scalarList alphasiOwn(alphas.size());
    scalarList alphasiNei(alphas.size());
    scalarList rhosiOwn(alphas.size());
    scalarList rhosiNei(alphas.size());

    scalarList alphaPhisi(alphas.size());
    scalarList alphaRhoPhisi(alphas.size());

    forAll(UOwn, facei)
    {
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
    const volScalarField& alpha1,
    const volScalarField& rho1,
    const volScalarField& rho2,
    const volVectorField& U,
    const volScalarField& e,
    const volScalarField& p,
    const volScalarField& c,
    surfaceScalarField& phi,
    surfaceScalarField& alphaPhi1,
    surfaceScalarField& alphaRhoPhi1,
    surfaceScalarField& alphaRhoPhi2,
    surfaceScalarField& rhoPhi,
    surfaceVectorField& rhoUPhi,
    surfaceScalarField& rhoEPhi
)
{
    createSavedFields();

    // Interpolate fields
    autoPtr<MUSCLReconstructionScheme<scalar>> alpha1Limiter
    (
        MUSCLReconstructionScheme<scalar>::New(alpha1, "alpha", alpha1.group())
    );
    autoPtr<MUSCLReconstructionScheme<scalar>> rho1Limiter
    (
        MUSCLReconstructionScheme<scalar>::New(rho1, "rho", rho1.group())
    );
    autoPtr<MUSCLReconstructionScheme<scalar>> rho2Limiter
    (
        MUSCLReconstructionScheme<scalar>::New(rho2, "rho", rho2.group())
    );
    autoPtr<MUSCLReconstructionScheme<vector>> ULimiter
    (
        MUSCLReconstructionScheme<vector>::New(U, "U")
    );
    autoPtr<MUSCLReconstructionScheme<scalar>> eLimiter
    (
        MUSCLReconstructionScheme<scalar>::New(e, "e")
    );
    autoPtr<MUSCLReconstructionScheme<scalar>> pLimiter
    (
        MUSCLReconstructionScheme<scalar>::New(p, "p")
    );
    autoPtr<MUSCLReconstructionScheme<scalar>> cLimiter
    (
        MUSCLReconstructionScheme<scalar>::New(c, "speedOfSound")
    );

    tmp<surfaceScalarField> talpha1Own(alpha1Limiter->interpolateOwn());
    tmp<surfaceScalarField> talpha1Nei(alpha1Limiter->interpolateNei());
    const surfaceScalarField& alpha1Own = talpha1Own();
    const surfaceScalarField& alpha1Nei = talpha1Nei();

    tmp<surfaceScalarField> talpha2Own(1.0 - alpha1Own);
    tmp<surfaceScalarField> talpha2Nei(1.0 - alpha1Nei);
    const surfaceScalarField& alpha2Own = talpha2Own();
    const surfaceScalarField& alpha2Nei = talpha2Nei();

    tmp<surfaceScalarField> trho1Own(rho1Limiter->interpolateOwn());
    tmp<surfaceScalarField> trho1Nei(rho1Limiter->interpolateNei());
    const surfaceScalarField& rho1Own = trho1Own();
    const surfaceScalarField& rho1Nei = trho1Nei();

    tmp<surfaceScalarField> trho2Own(rho2Limiter->interpolateOwn());
    tmp<surfaceScalarField> trho2Nei(rho2Limiter->interpolateNei());
    const surfaceScalarField& rho2Own = trho2Own();
    const surfaceScalarField& rho2Nei = trho2Nei();

    tmp<surfaceVectorField> tUOwn(ULimiter->interpolateOwn());
    tmp<surfaceVectorField> tUNei(ULimiter->interpolateNei());
    const surfaceVectorField& UOwn = tUOwn();
    const surfaceVectorField& UNei = tUNei();

    tmp<surfaceScalarField> teOwn(eLimiter->interpolateOwn());
    tmp<surfaceScalarField> teNei(eLimiter->interpolateNei());
    const surfaceScalarField& eOwn = teOwn();
    const surfaceScalarField& eNei = teNei();

    tmp<surfaceScalarField> tpOwn(pLimiter->interpolateOwn());
    tmp<surfaceScalarField> tpNei(pLimiter->interpolateNei());
    const surfaceScalarField& pOwn = tpOwn();
    const surfaceScalarField& pNei = tpNei();

    tmp<surfaceScalarField> tcOwn(cLimiter->interpolateOwn());
    tmp<surfaceScalarField> tcNei(cLimiter->interpolateNei());
    const surfaceScalarField& cOwn = tcOwn();
    const surfaceScalarField& cNei = tcNei();

    surfaceScalarField rhoOwn(alpha1Own*rho1Own + alpha2Own*rho2Own);
    surfaceScalarField rhoNei(alpha1Nei*rho1Nei + alpha2Nei*rho2Nei);

    preUpdate(p);

    forAll(UOwn, facei)
    {
        scalarList alphaPhisi(2);
        scalarList alphaRhoPhisi(2);
        calculateFluxes
        (
            {alpha1Own[facei], alpha2Own[facei]},
            {alpha1Nei[facei], alpha2Nei[facei]},
            {rho1Own[facei], rho2Own[facei]},
            {rho1Nei[facei], rho2Nei[facei]},
            rhoOwn[facei], rhoNei[facei],
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

        alphaPhi1[facei] = alphaPhisi[0];
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
                    alpha1Own.boundaryField()[patchi][facei],
                    alpha2Own.boundaryField()[patchi][facei]
                },
                {
                    alpha1Nei.boundaryField()[patchi][facei],
                    alpha2Nei.boundaryField()[patchi][facei]
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
                mesh_.Sf().boundaryField()[patchi][facei],
                phi.boundaryFieldRef()[patchi][facei],
                alphaPhisi,
                alphaRhoPhisi,
                rhoUPhi.boundaryFieldRef()[patchi][facei],
                rhoEPhi.boundaryFieldRef()[patchi][facei],
                facei, patchi
            );

            alphaPhi1.boundaryFieldRef()[patchi][facei] = alphaPhisi[0];
            alphaRhoPhi1.boundaryFieldRef()[patchi][facei] =
                alphaRhoPhisi[0];
            alphaRhoPhi2.boundaryFieldRef()[patchi][facei] =
                alphaRhoPhisi[1];

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
    autoPtr<MUSCLReconstructionScheme<scalar>> rhoLimiter
    (
        MUSCLReconstructionScheme<scalar>::New(rho, "rho")
    );
    tmp<surfaceScalarField> trhoOwn(rhoLimiter->interpolateOwn());
    tmp<surfaceScalarField> trhoNei(rhoLimiter->interpolateNei());
    surfaceScalarField& rhoOwn = trhoOwn.ref();
    surfaceScalarField& rhoNei = trhoNei.ref();

    // Interpolate fields
    autoPtr<MUSCLReconstructionScheme<vector>> ULimiter
    (
        MUSCLReconstructionScheme<vector>::New(U, "U")
    );
    autoPtr<MUSCLReconstructionScheme<scalar>> eLimiter
    (
        MUSCLReconstructionScheme<scalar>::New(e, "e")
    );
    autoPtr<MUSCLReconstructionScheme<scalar>> pLimiter
    (
        MUSCLReconstructionScheme<scalar>::New(p, "p")
    );

    tmp<surfaceVectorField> tUOwn(ULimiter->interpolateOwn());
    tmp<surfaceVectorField> tUNei(ULimiter->interpolateNei());
    const surfaceVectorField& UOwn = tUOwn();
    const surfaceVectorField& UNei = tUNei();

    tmp<surfaceScalarField> teOwn(eLimiter->interpolateOwn());
    tmp<surfaceScalarField> teNei(eLimiter->interpolateNei());
    const surfaceScalarField& eOwn = teOwn();
    const surfaceScalarField& eNei = teNei();

    tmp<surfaceScalarField> tpOwn(pLimiter->interpolateOwn());
    tmp<surfaceScalarField> tpNei(pLimiter->interpolateNei());
    const surfaceScalarField& pOwn = tpOwn();
    const surfaceScalarField& pNei = tpNei();

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
            rhoOwn[facei], rhoNei[facei],
            UOwn[facei], UNei[facei],
            eOwn[facei], eNei[facei],
            pOwn[facei], pNei[facei],
            mesh_.Sf()[facei],
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
                    rhoOwn.boundaryField()[patchi][facei],
                    rhoNei.boundaryField()[patchi][facei],
                    UOwn.boundaryField()[patchi][facei],
                    UNei.boundaryField()[patchi][facei],
                    eOwn.boundaryField()[patchi][facei],
                    eNei.boundaryField()[patchi][facei],
                    pOwn.boundaryField()[patchi][facei],
                    pNei.boundaryField()[patchi][facei],
                    mesh_.Sf().boundaryField()[patchi][facei],
                    facei, patchi
                );
        }
    }
    return tmpPhi;
}

// ************************************************************************* //
