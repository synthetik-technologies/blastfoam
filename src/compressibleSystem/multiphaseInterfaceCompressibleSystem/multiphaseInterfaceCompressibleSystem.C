/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "multiphaseInterfaceCompressibleSystem.H"
// #include "alphaContactAngleFvPatchScalarField.H"
#include "extendedNLevelGlobalCellToCellStencils.H"
#include "MULES.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiphaseInterfaceCompressibleSystem, 0);
    addToRunTimeSelectionTable
    (
        compressibleSystem,
        multiphaseInterfaceCompressibleSystem,
        multiphase
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseInterfaceCompressibleSystem::multiphaseInterfaceCompressibleSystem
(
    const fvMesh& mesh
)
:
    multiphaseCompressibleSystem(mesh),
    interfaceSystem(U_, *this)
{
    ITstream is(this->lookup("sigmas"));
    token t(is);
    if (!t.isPunctuation() || t.pToken() != token::BEGIN_LIST)
    {
        FatalIOErrorInFunction(is)
            << "Expected " << token::BEGIN_LIST << " but found " << t.info() << endl
            << abort(FatalIOError);
    }
    while (is.good())
    {
        is  >> t;
        if (t.isPunctuation() && t.pToken() == token::END_LIST)
        {
            break;
        }
        is.putBack(t);
        interfacePair pair(is);
        autoPtr<surfaceTensionModel> sfPtr =
            surfaceTensionModel::New
            (
                dictionary(is),
                mesh
            );
        sfPtr->rename
        (
            IOobject::groupName
            (
                surfaceTensionModel::typeName,
                pair.name()
            )
        );
        stModels_.insert(pair,sfPtr.ptr());
    }
}


Foam::multiphaseInterfaceCompressibleSystem::multiphaseInterfaceCompressibleSystem
(
    const fvMesh& mesh,
    const bool
)
:
    multiphaseCompressibleSystem(mesh),
    interfaceSystem(U_, *this)
{
    ITstream is(this->lookup("sigmas"));
    token t(is);
    if (!t.isPunctuation() || t.pToken() != token::BEGIN_LIST)
    {
        FatalIOErrorInFunction(is)
            << "Expected " << token::BEGIN_LIST << " but found " << t.info() << endl
            << abort(FatalIOError);
    }
    while (is.good())
    {
        is  >> t;
        if (t.isPunctuation() && t.pToken() == token::END_LIST)
        {
            break;
        }
        is.putBack(t);
        interfacePair pair(is);
        autoPtr<surfaceTensionModel> sfPtr =
            surfaceTensionModel::New
            (
                dictionary(is),
                mesh
            );
        sfPtr->rename
        (
            IOobject::groupName
            (
                surfaceTensionModel::typeName,
                pair.name()
            )
        );
        stModels_.insert(pair,sfPtr.ptr());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiphaseInterfaceCompressibleSystem::~multiphaseInterfaceCompressibleSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multiphaseInterfaceCompressibleSystem::update()
{
    decode();

    PtrList<surfaceScalarField> rhosOwn(rhos_.size());
    PtrList<surfaceScalarField> rhosNei(rhos_.size());
    PtrList<surfaceScalarField> alphasOwn(rhos_.size());
    PtrList<surfaceScalarField> alphasNei(rhos_.size());

    surfaceScalarField rhoOwn
    (
        surfaceScalarField::New
        (
            "rhoOwn",
            mesh(),
            dimensionedScalar(dimDensity, 0.0)
        )
    );
    surfaceScalarField rhoNei
    (
        surfaceScalarField::New
        (
            "rhoNei",
            mesh(),
            dimensionedScalar(dimDensity, 0.0)
        )
    );
    forAll(rhos_, phasei)
    {
        autoPtr<ReconstructionScheme<scalar>> alphaLimiter
        (
            ReconstructionScheme<scalar>::New
            (
                alphas_[phasei],
                "alpha",
                alphas_[phasei].group(),
                true
            )
        );
        alphasOwn.set(phasei, alphaLimiter->interpolateOwn());
        alphasNei.set(phasei, alphaLimiter->interpolateNei());

        autoPtr<ReconstructionScheme<scalar>> rhoLimiter
        (
            ReconstructionScheme<scalar>::New
            (
                rhos_[phasei],
                "rho",
                rhos_[phasei].group(),
                true
            )
        );
        rhosOwn.set(phasei, rhoLimiter->interpolateOwn());
        rhosNei.set(phasei, rhoLimiter->interpolateNei());
        fluxScheme::correctPhaseFields
        (
            alphas_[phasei],
            rhosOwn[phasei], rhosNei[phasei],
            thermo_.thermo(phasei).residualAlpha().value()
        );

        tmp<surfaceScalarField> talphaRhoOwn = surfaceScalarField::New
        (
            rhoLimiter->ownName(alphaRhos_[phasei].name()),
            alphasOwn[phasei]*rhosOwn[phasei]
        );
        tmp<surfaceScalarField> talphaRhoNei = surfaceScalarField::New
        (
            rhoLimiter->neiName(alphaRhos_[phasei].name()),
            alphasNei[phasei]*rhosNei[phasei]
        );

        if (mesh().cacheTemporaryObject(talphaRhoOwn().name()))
        {
            mesh().cacheTemporaryObject(talphaRhoOwn.ref());
            mesh().cacheTemporaryObject(talphaRhoNei.ref());
        }

        rhoOwn += talphaRhoOwn;
        rhoNei += talphaRhoNei;
    }
    if (mesh().cacheTemporaryObject(rhoOwn.name()))
    {
        mesh().cacheTemporaryObject(rhoOwn);
        mesh().cacheTemporaryObject(rhoNei);
    }

    fluxScheme_->update
    (
        rhoOwn,
        rhoNei,
        U_,
        e_,
        p_,
        speedOfSound()(),
        phi_,
        rhoPhi_,
        rhoUPhi_,
        rhoEPhi_
    );


    // Limit alpha flux
    volScalarField divPhi(fvc::div(phi_));
    surfaceScalarField phi(phi_);
    this->storeAndBlendDelta(phi);

    UPtrList<const volScalarField> alphas(alphas_.size());
    forAll(alphas_, phasei)
    {
        alphas.set(phasei, &alphas_[phasei]);
        volScalarField alphaOld(alphas_[phasei]);
        this->blendOld(alphaOld);
        alphaOld.storeOldTime();

        surfaceScalarField& alphaPhi = alphaPhis_[phasei];
        alphaPhi =
            fluxScheme_->flux(alphasOwn[phasei], alphasNei[phasei], phi_);
        this->storeAndBlendDelta(alphaPhi);

        MULES::limit
        (
            1.0/mesh().time().deltaT().value(),
            geometricOneField(),
            alphaOld,
            phi,
            alphaPhi,
            zeroField(),
            zeroField(),//(-divPhi*alphas_[phasei])(),
            oneField(),
            zeroField(),
            false
        );
        alphaPhi = this->calcAndStoreDelta(alphaPhi);
    }
    MULES::limitSum(alphas, alphaPhis_, phi_);

    // Update phase mass fluxes
    forAll(alphas_, phasei)
    {
        alphaRhoPhis_[phasei] =
            fluxScheme_->flux
            (
                rhosOwn[phasei],
                rhosNei[phasei],
                alphaPhis_[phasei]
            );
    }

    PtrList<surfaceScalarField> alphaRhoPhiUDs(alphaRhoPhis_.size());
    forAll(alphaRhoPhis_, phasei)
    {
        alphaRhoPhiUDs.set
        (
            phasei,
            upwind<scalar>(mesh(), alphaPhis_[phasei]).flux(rhos_[phasei])
        );

        alphaRhoPhis_[phasei] -= alphaRhoPhiUDs[phasei];
    }

    {
        UPtrList<scalarField> alphaRhoPhisInternal(alphaRhoPhis_.size());

        forAll(alphaRhoPhisInternal, phasei)
        {
            alphaRhoPhisInternal.set(phasei, &alphaRhoPhis_[phasei]);
        }

        MULES::limitSum(alphaRhoPhisInternal);
    }

    const surfaceScalarField::Boundary& phibf = phi_.boundaryField();
    forAll(phibf, patchi)
    {
        if (phibf[patchi].coupled())
        {
            UPtrList<scalarField> alphaRhoPhisPatch(alphaRhoPhis_.size());

            forAll(alphaRhoPhisPatch, phasei)
            {
                alphaRhoPhisPatch.set
                (
                    phasei,
                    &alphaRhoPhis_[phasei].boundaryFieldRef()[patchi]
                );
            }

            MULES::limitSum(alphaRhoPhisPatch);
        }
    }

    forAll(alphaRhoPhis_, phasei)
    {
        alphaRhoPhis_[phasei] += alphaRhoPhiUDs[phasei];
    }

    thermo_.update();
}


// void Foam::multiphaseInterfaceCompressibleSystem::decode()
// {
//     rho_ == Zero;
//     const extendedNLevelCFCCellToCellStencil& stencil
//     (
//         extendedNLevelCFCCellToCellStencil::New(mesh(), 3)
//     );
//     List<List<scalar>> alphaNei(rho_.size());
//     List<List<scalar>> rhoNei(rho_.size());
//
//     forAll(alphas_, phasei)
//     {
//         volScalarField& alpha = alphas_[phasei];
//         volScalarField& rho = rhos_[phasei];
//         volScalarField& alphaRho = alphaRhos_[phasei];
//
//         alpha.maxMin(0.0, 1.0);
//         alphaRho.max(0);
//
//
//         // Only update cells that have a valid volume fraction
//         // other cell densities are handled by transport of density
//         const scalar rAlpha = thermo_.thermo(phasei).residualAlpha().value();
//         if (usesCompression(alpha.name()))
//         {
//             forAll(alpha, celli)
//             {
//                 if (alpha[celli] > rAlpha)
//                 {
//                     rho[celli] = alphaRho[celli]/alpha[celli];
//                 }
//             }
//
//             stencil.collectData(alpha, alphaNei);
//             stencil.collectData(rho, rhoNei);
//
//             forAll(alpha, celli)
//             {
//                 const scalar alphai = alpha[celli];
//                 if (alphai < 0.5)
//                 {
//                     scalar sumRhoW = 0.0;
//                     scalar sumW = 0.0;
//                     forAll(rhoNei[celli], cj)
//                     {
//                         if (alphaNei[celli][cj] > 0.5)
//                         {
//                             scalar w = 1.0/max(1.0 - alphaNei[celli][cj], rAlpha);
//                             sumRhoW += rhoNei[celli][cj]*w;
//                             sumW += w;
//                         }
//                     }
//                     rho[celli] = sumRhoW/max(sumW, rAlpha);
//                 }
//             }
//         }
//         else
//         {
//             rho.ref() = alphaRho()/max(alpha(), rAlpha);
//         }
//         rho.correctBoundaryConditions();
//         alphaRho.correctBoundaryConditions();
//         alphaRho.boundaryFieldRef() = alpha.boundaryField()*rho.boundaryField();
//         rho_ += alphaRho;
//     }
//
//     compressibleBlastSystem::decode();
// }

// ************************************************************************* //
