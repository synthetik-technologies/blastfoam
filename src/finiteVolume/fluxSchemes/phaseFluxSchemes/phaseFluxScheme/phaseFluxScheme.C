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

#include "phaseFluxScheme.H"
#include "MUSCLReconstructionScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseFluxScheme, 0);
    defineRunTimeSelectionTable(phaseFluxScheme, dictionary);
    defineRunTimeSelectionTable(phaseFluxScheme, solid);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseFluxScheme::phaseFluxScheme(const fvMesh& mesh, const word& name)
:
    fluxSchemeBase(mesh, name),
    dict_(mesh.schemesDict().subDict("fluxSchemes").subDict(name))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseFluxScheme::~phaseFluxScheme()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::phaseFluxScheme::interpolate
(
    const volScalarField& f,
    const word& fName
) const
{
    return interpolateField(f, fName);
}


Foam::tmp<Foam::surfaceVectorField> Foam::phaseFluxScheme::interpolate
(
    const volVectorField& f,
    const word& fName
) const
{
    return interpolateField(f, fName);
}


Foam::tmp<Foam::surfaceSymmTensorField> Foam::phaseFluxScheme::interpolate
(
    const volSymmTensorField& f,
    const word& fName
) const
{
    return interpolateField(f, fName);
}


Foam::tmp<Foam::surfaceSphericalTensorField> Foam::phaseFluxScheme::interpolate
(
    const volSphericalTensorField& f,
    const word& fName
) const
{
    return interpolateField(f, fName);
}


Foam::tmp<Foam::surfaceTensorField> Foam::phaseFluxScheme::interpolate
(
    const volTensorField& f,
    const word& fName
) const
{
    return interpolateField(f, fName);
}


void Foam::phaseFluxScheme::clear()
{
    Uf_.clear();
    pf_.clear();
    alphaf_.clear();
}

void Foam::phaseFluxScheme::createSavedFields()
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
                IOobject::groupName("fluxScheme::Uf", this->group()),
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedVector("0", dimVelocity, Zero)
        )
    );
    pf_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                IOobject::groupName("fluxScheme::pf", this->group()),
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimPressure, Zero)
        )
    );
    alphaf_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                IOobject::groupName("phaseFluxScheme::alphaf", this->group()),
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimless, Zero)
        )
    );
}

Foam::tmp<Foam::surfaceVectorField> Foam::phaseFluxScheme::Uf() const
{
    if (Uf_.valid())
    {
        return Uf_();
    }
    FatalErrorInFunction
        << IOobject::groupName("phaseFluxScheme::Uf", this->group())
        << " has not been set." << nl
        << abort(FatalError);

    return Uf_();
}


Foam::tmp<Foam::surfaceScalarField> Foam::phaseFluxScheme::pf() const
{
    if (pf_.valid())
    {
        return pf_();
    }
    FatalErrorInFunction
        << IOobject::groupName("phaseFluxScheme::pf", this->group())
        << " has not been set." << nl
        << abort(FatalError);

    return pf_;
}


Foam::tmp<Foam::surfaceScalarField> Foam::phaseFluxScheme::alphaf() const
{
    if (alphaf_.valid())
    {
        return alphaf_();
    }
    FatalErrorInFunction
        << IOobject::groupName("phaseFluxScheme::alphaf", this->group())
        << " has not been set." << nl
        << abort(FatalError);

    return alphaf_;
}


void Foam::phaseFluxScheme::update
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& e,
    const volScalarField& p,
    const volScalarField& c,
    surfaceScalarField& phi,
    surfaceScalarField& alphaPhi,
    surfaceScalarField& alphaRhoPhi,
    surfaceVectorField& alphaRhoUPhi,
    surfaceScalarField& alphaRhoEPhi
)
{
    createSavedFields();

    autoPtr<MUSCLReconstructionScheme<scalar>> alphaLimiter
    (
        MUSCLReconstructionScheme<scalar>::New(alpha, "alpha")
    );
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

    tmp<surfaceScalarField> talphaOwn(alphaLimiter->interpolateOwn());
    tmp<surfaceScalarField> talphaNei(alphaLimiter->interpolateNei());
    const surfaceScalarField& alphaOwn = talphaOwn();
    const surfaceScalarField& alphaNei = talphaNei();

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
        if (alphaOwn[facei] < 1e-10 && alphaNei[facei] < 1e-10)
        {
            phi[facei] = 0.0;
            alphaPhi[facei] = 0.0;
            alphaRhoPhi[facei] = 0.0;
            alphaRhoUPhi[facei] = Zero;
            alphaRhoEPhi[facei] = 0.0;
        }
        else
        {
            calculateFluxes
            (
                alphaOwn[facei], alphaNei[facei],
                rhoOwn[facei], rhoNei[facei],
                UOwn[facei], UNei[facei],
                eOwn[facei], eNei[facei],
                pOwn[facei], pNei[facei],
                cOwn[facei], cNei[facei],
                mesh_.Sf()[facei],
                phi[facei],
                alphaPhi[facei],
                alphaRhoPhi[facei],
                alphaRhoUPhi[facei],
                alphaRhoEPhi[facei],
                facei
            );
        }
    }

    forAll(U.boundaryField(), patchi)
    {
        forAll(U.boundaryField()[patchi], facei)
        {
            if
            (
                alphaOwn.boundaryField()[patchi][facei] < 1e-10
             && alphaNei.boundaryField()[patchi][facei] < 1e-10
            )
            {
                phi.boundaryFieldRef()[patchi][facei] = 0.0;
                alphaPhi.boundaryFieldRef()[patchi][facei] = 0.0;
                alphaRhoPhi.boundaryFieldRef()[patchi][facei] = 0.0;
                alphaRhoUPhi.boundaryFieldRef()[patchi][facei] = Zero;
                alphaRhoEPhi.boundaryFieldRef()[patchi][facei] = 0.0;
            }
            else
            {
                calculateFluxes
                (
                    alphaOwn.boundaryField()[patchi][facei],
                    alphaNei.boundaryField()[patchi][facei],
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
                    alphaPhi.boundaryFieldRef()[patchi][facei],
                    alphaRhoPhi.boundaryFieldRef()[patchi][facei],
                    alphaRhoUPhi.boundaryFieldRef()[patchi][facei],
                    alphaRhoEPhi.boundaryFieldRef()[patchi][facei],
                    facei, patchi
                );
            }
        }
    }
    postUpdate();
}

void Foam::phaseFluxScheme::update
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& e,
    const volScalarField& p,
    const volScalarField& c,
    surfaceScalarField& phi,
    surfaceScalarField& alphaRhoPhi,
    surfaceVectorField& alphaRhoUPhi,
    surfaceScalarField& alphaRhoEPhi
)
{
    createSavedFields();
    const word phaseName(U.group());

    autoPtr<MUSCLReconstructionScheme<scalar>> alphaLimiter
    (
        MUSCLReconstructionScheme<scalar>::New(alpha, "alpha", phaseName)
    );
    autoPtr<MUSCLReconstructionScheme<scalar>> rhoLimiter
    (
        MUSCLReconstructionScheme<scalar>::New(rho, "rho", phaseName)
    );
    autoPtr<MUSCLReconstructionScheme<vector>> ULimiter
    (
        MUSCLReconstructionScheme<vector>::New(U, "U", phaseName)
    );
    autoPtr<MUSCLReconstructionScheme<scalar>> eLimiter
    (
        MUSCLReconstructionScheme<scalar>::New(e, "e", phaseName)
    );
    autoPtr<MUSCLReconstructionScheme<scalar>> pLimiter
    (
        MUSCLReconstructionScheme<scalar>::New(p, "p", phaseName)
    );
    autoPtr<MUSCLReconstructionScheme<scalar>> cLimiter
    (
        MUSCLReconstructionScheme<scalar>::New(c, "speedOfSound", phaseName)
    );

    tmp<surfaceScalarField> talphaOwn(alphaLimiter->interpolateOwn());
    tmp<surfaceScalarField> talphaNei(alphaLimiter->interpolateNei());
    const surfaceScalarField& alphaOwn = talphaOwn();
    const surfaceScalarField& alphaNei = talphaNei();

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
        scalar alphaPhi;
        calculateFluxes
        (
            alphaOwn[facei], alphaNei[facei],
            rhoOwn[facei], rhoNei[facei],
            UOwn[facei], UNei[facei],
            eOwn[facei], eNei[facei],
            pOwn[facei], pNei[facei],
            cOwn[facei], cNei[facei],
            mesh_.Sf()[facei],
            phi[facei],
            alphaPhi,
            alphaRhoPhi[facei],
            alphaRhoUPhi[facei],
            alphaRhoEPhi[facei],
            facei
        );
    }

    forAll(U.boundaryField(), patchi)
    {
        forAll(U.boundaryField()[patchi], facei)
        {
            if
            (
                alphaOwn.boundaryField()[patchi][facei] < 1e-10
             && alphaNei.boundaryField()[patchi][facei] < 1e-10
            )
            {
                continue;
            }
            scalar alphaPhi;
            calculateFluxes
            (
                alphaOwn.boundaryField()[patchi][facei],
                alphaNei.boundaryField()[patchi][facei],
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
                alphaPhi,
                alphaRhoPhi.boundaryFieldRef()[patchi][facei],
                alphaRhoUPhi.boundaryFieldRef()[patchi][facei],
                alphaRhoEPhi.boundaryFieldRef()[patchi][facei],
                facei, patchi
            );
        }
    }
    postUpdate();
}


void Foam::phaseFluxScheme::update
(
    const UPtrList<volScalarField>& alphas,
    const UPtrList<volScalarField>& rhos,
    const volVectorField& U,
    const volScalarField& e,
    const volScalarField& p,
    const volScalarField& c,
    surfaceScalarField& phi,
    surfaceScalarField& alphaPhi,
    surfaceScalarField& alphaRhoPhi,
    PtrList<surfaceScalarField>& alphaPhis,
    PtrList<surfaceScalarField>& alphaRhoPhis,
    surfaceVectorField& alphaRhoUPhi,
    surfaceScalarField& alphaRhoEPhi
)
{
    createSavedFields();
    const word phaseName(U.group());

    // Interpolate fields
    PtrList<surfaceScalarField> alphasOwn(alphas.size());
    PtrList<surfaceScalarField> alphasNei(alphas.size());

    PtrList<surfaceScalarField> rhosOwn(alphas.size());
    PtrList<surfaceScalarField> rhosNei(alphas.size());

    surfaceScalarField alphaOwn
    (
        IOobject
        (
            IOobject::groupName("alphaOwn", phaseName),
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("0", dimless, 0.0)
    );
    surfaceScalarField alphaNei
    (
        IOobject
        (
            IOobject::groupName("alphaNei", phaseName),
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("0", dimless, 0.0)
    );

    surfaceScalarField rhoOwn
    (
        IOobject
        (
            IOobject::groupName("rhoNei", phaseName),
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("0", dimDensity, 0.0)
    );
    surfaceScalarField rhoNei
    (
        IOobject
        (
            IOobject::groupName("rhoNei", phaseName),
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("0", dimDensity, 0.0)
    );

    forAll(alphas, phasei)
    {
        const word phaseNamei(alphas[phasei].group());
        autoPtr<MUSCLReconstructionScheme<scalar>> alphaLimiter
        (
            MUSCLReconstructionScheme<scalar>::New
            (
                alphas[phasei],
                "alpha",
                phaseNamei
            )
        );
        autoPtr<MUSCLReconstructionScheme<scalar>> rhoLimiter
        (
            MUSCLReconstructionScheme<scalar>::New(rhos[phasei], "rho", phaseNamei)
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
        alphaOwn += alphasOwn[phasei];
        alphaNei += alphasNei[phasei];
        rhoOwn += alphasOwn[phasei]*rhosOwn[phasei];
        rhoNei += alphasNei[phasei]*rhosNei[phasei];
    }
    rhoOwn /= max(alphaOwn, 1e-10);
    rhoNei /= max(alphaNei, 1e-10);

    autoPtr<MUSCLReconstructionScheme<vector>> ULimiter
    (
        MUSCLReconstructionScheme<vector>::New(U, "U", phaseName)
    );
    autoPtr<MUSCLReconstructionScheme<scalar>> eLimiter
    (
        MUSCLReconstructionScheme<scalar>::New(e, "e", phaseName)
    );
    autoPtr<MUSCLReconstructionScheme<scalar>> pLimiter
    (
        MUSCLReconstructionScheme<scalar>::New(p, "p", phaseName)
    );
    autoPtr<MUSCLReconstructionScheme<scalar>> cLimiter
    (
        MUSCLReconstructionScheme<scalar>::New(c, "speedOfSound", phaseName)
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
    forAll(UOwn, facei)
    {
        if (alphaOwn[facei] < small && alphaNei[facei] < small)
        {
            continue;
        }
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
            alphaOwn[facei], alphaNei[facei],
            rhoOwn[facei], rhoNei[facei],
            alphasiOwn, alphasiNei,
            rhosiOwn, rhosiNei,
            UOwn[facei], UNei[facei],
            eOwn[facei], eNei[facei],
            pOwn[facei], pNei[facei],
            cOwn[facei], cNei[facei],
            mesh_.Sf()[facei],
            phi[facei],
            alphaPhisi,
            alphaRhoPhisi,
            alphaRhoUPhi[facei],
            alphaRhoEPhi[facei],
            facei
        );

        alphaRhoPhi[facei] = 0.0;
        alphaPhi[facei] = 0.0;
        forAll(alphas, phasei)
        {
            alphaPhis[phasei][facei] = alphaPhisi[phasei];
            alphaRhoPhis[phasei][facei] = alphaRhoPhisi[phasei];
            alphaRhoPhi[facei] += alphaRhoPhisi[phasei];
            alphaPhi[facei] += alphaPhisi[phasei];
        }
    }

    forAll(U.boundaryField(), patchi)
    {
        forAll(U.boundaryField()[patchi], facei)
        {
            if
            (
                alphaOwn.boundaryField()[patchi][facei] < small
             && alphaNei.boundaryField()[patchi][facei] < small)
            {
                continue;
            }
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
                alphaOwn.boundaryField()[patchi][facei],
                alphaNei.boundaryField()[patchi][facei],
                rhoOwn.boundaryField()[patchi][facei],
                rhoNei.boundaryField()[patchi][facei],
                alphasiOwn,
                alphasiNei,
                rhosiOwn,
                rhosiNei,
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
                alphaRhoUPhi.boundaryFieldRef()[patchi][facei],
                alphaRhoEPhi.boundaryFieldRef()[patchi][facei],
                facei, patchi
            );

            alphaPhi.boundaryFieldRef()[patchi][facei] = 0.0;
            alphaRhoPhi.boundaryFieldRef()[patchi][facei] = 0.0;
            forAll(alphas, phasei)
            {
                alphaPhis[phasei].boundaryFieldRef()[patchi][facei] =
                    alphaPhisi[phasei];
                alphaRhoPhis[phasei].boundaryFieldRef()[patchi][facei] =
                    alphaRhoPhisi[phasei];
                alphaPhi.boundaryFieldRef()[patchi][facei] +=
                    alphaPhisi[phasei];
                alphaRhoPhi.boundaryFieldRef()[patchi][facei] +=
                    alphaRhoPhisi[phasei];
            }
        }
    }
    postUpdate();
}

// ************************************************************************* //
