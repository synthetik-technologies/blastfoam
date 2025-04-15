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
#include "dynamicMultiAlphaContactAngleFvPatchScalarField.H"
#include "extendedNLevelGlobalCellToCellStencils.H"
#include "unitConversion.H"
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

// * * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * //

// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::multiphaseInterfaceCompressibleSystem::correctContactAngle
(
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    surfaceVectorField::Boundary& nHatfv
) const
{
    typedef dynamicMultiAlphaContactAngleFvPatchScalarField multiAlphaContact;
    const volScalarField::Boundary& a1bf = alpha1.boundaryField();
    const volScalarField::Boundary& a2bf = alpha2.boundaryField();

    const fvBoundaryMesh& boundary = mesh().boundary();

    forAll(boundary, patchi)
    {
        if
        (
            isA<multiAlphaContact>(a1bf[patchi])
         || isA<multiAlphaContact>(a2bf[patchi])
        )
        {
            if
            (
                isA<multiAlphaContact>(a1bf[patchi])
             && isA<multiAlphaContact>(a2bf[patchi])
            )
            {
                FatalErrorInFunction
                    << "alphaContactAngle boundary condition "
                       "specified on patch " << boundary[patchi].name()
                    << " for both " << alpha1.name() << " and " << alpha2.name()
                    << nl << "which may be inconsistent."
                    << exit(FatalError);
            }

            const fvPatch& patch = boundary[patchi];

            const multiAlphaContact& acap =
                isA<multiAlphaContact>(a1bf[patchi])
              ? refCast<const multiAlphaContact>(a1bf[patchi])
              : refCast<const multiAlphaContact>(a2bf[patchi])
              ;

            vectorField& nHatPatch = nHatfv[patchi];
            const vectorField nf(patch.nf());

            multiAlphaContact::thetaPropsTable::
                const_iterator tp =
                acap.thetaProps().find(interfacePair(alpha1, alpha2));

            if (tp == acap.thetaProps().end())
            {
                FatalErrorInFunction
                    << "Cannot find interface " << interfacePair(alpha1, alpha2)
                    << "\n    in table of theta properties for patch "
                    << acap.patch().name()
                    << exit(FatalError);
            }

            const bool matched = (tp.key().first() == alpha1.name());

            const scalar theta0 = degToRad(tp().theta0(matched));

            scalarField theta(boundary[patchi].size(), theta0);

            const scalar uTheta = tp().uTheta();

            // Calculate the dynamic contact angle if required
            if (uTheta > small)
            {
                const scalar thetaA = degToRad(tp().thetaA(matched));
                const scalar thetaR = degToRad(tp().thetaR(matched));

                // Calculated the component of the velocity parallel to the wall
                vectorField Uwall
                (
                    U_.boundaryField()[patchi].patchInternalField()
                  - U_.boundaryField()[patchi]
                );
                Uwall -= (nf & Uwall)*nf;

                // Find the direction of the interface parallel to the wall
                vectorField nWall(nHatPatch - (nf & nHatPatch)*nf);

                // Normalise nWall
                nWall /= (mag(nWall) + small);

                // Calculate Uwall resolved normal to the interface parallel to
                // the interface
                const scalarField uwall(nWall & Uwall);

                theta += (thetaA - thetaR)*tanh(uwall/uTheta);
            }


            // Reset nHatPatch to correspond to the contact angle

            forAll(nHatPatch, facei)
            {
                const scalar a12 = nHatPatch[facei] & nf[facei];
                const scalar b1 = cos(theta[facei]);
                const scalar b2 = cos(acos(a12) - theta[facei]);
                const scalar det(1.0 - a12*a12);
                const scalar a((b1 - a12*b2)/det);
                const scalar b((b2 - a12*b1)/det);

                nHatPatch[facei] = a*nf[facei] + b*nHatPatch[facei];

                nHatPatch[facei] /= (mag(nHatPatch[facei]) + deltaN_.value());
            }
        }
        else if (isA<alphaContactAngleFvPatchScalarField>(a1bf[patchi]))
        {
            FatalErrorInFunction
                << " When using more that 2 phases the "
                << multiAlphaContact::typeName
                << " boundary condition should be used in place of "
                << a1bf[patchi].type() << endl
                << abort(FatalError);
        }
        else if (isA<alphaContactAngleFvPatchScalarField>(a2bf[patchi]))
        {
            FatalErrorInFunction
                << " When using more that 2 phases the "
                << multiAlphaContact::typeName
                << " boundary condition should be used in place of "
                << a2bf[patchi].type() << endl
                << abort(FatalError);
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInterfaceCompressibleSystem::K
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    surfaceVectorField gradAlphaf
    (
        fvc::interpolate(alpha2)*fvc::interpolate(fvc::grad(alpha1))
      - fvc::interpolate(alpha1)*fvc::interpolate(fvc::grad(alpha2))
    );

    // Face unit interface normal
    tmp<surfaceVectorField> tnHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));

    correctContactAngle(alpha1, alpha2, tnHatfv.ref().boundaryFieldRef());

    // Simple expression for curvature
    return -fvc::div(tnHatfv & mesh().Sf());
}


void Foam::multiphaseInterfaceCompressibleSystem::addSources
(
    volVectorField::Internal& rhoUSource,
    volScalarField::Internal& rhoESource
) const
{

    multiphaseCompressibleSystem::addSources(rhoUSource, rhoESource);

    tmp<surfaceScalarField> tstf
    (
        surfaceScalarField::New
        (
            "surfaceTensionForce",
            mesh(),
            dimensionedScalar(dimensionSet(1, -2, -2, 0, 0), 0.0)
        )
    );
    surfaceScalarField& stf = tstf.ref();

    forAll(alphas_, phasei)
    {
        const volScalarField& alpha1 = alphas_[phasei];
        forAll(alphas_, phasej)
        {
            if (phasei == phasej)
            {
                continue;
            }
            const volScalarField& alpha2 = alphas_[phasej];
            typename surfaceTensionTable::const_iterator iter =
                stModels_.find(interfacePair(alpha1, alpha2));
            if (iter != stModels_.cend())
            {
                stf +=
                    fvc::interpolate(iter()->sigma()*K(alpha1, alpha2))
                   *(
                        fvc::interpolate(alpha2)*fvc::snGrad(alpha1)
                     - fvc::interpolate(alpha1)*fvc::snGrad(alpha2)
                    );
            }
        }
    }

    tmp<volVectorField> stF(fvc::reconstruct(tstf*mesh().magSf()));

    rhoUSource -= stF();
    rhoESource -= stF() & U_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseInterfaceCompressibleSystem::multiphaseInterfaceCompressibleSystem
(
    const fvMesh& mesh,
    const bool initialize
)
:
    multiphaseCompressibleSystem(mesh, initialize),
    deltaN_
    (
        "deltaN",
        1e-8/cbrt(min(this->mesh().V()))
    )
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
        interfacePair key(is);
        if (whichPhase(key.first()) >= 0 && whichPhase(key.second()) >= 0)
        {
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
                    key.name()
                )
            );
            stModels_.insert(key, sfPtr.ptr());
        }
        else
        {
            WarningInFunction
                << "Unknown phase is pair " << key <<endl;
        }
    }

    forAll(alphas_, phasei)
    {
        const volScalarField& alpha1 = alphas_[phasei];
        for (label phasej = phasei+1; phasej < alphas_.size(); phasej++)
        {
            const volScalarField& alpha2 = alphas_[phasej];
            typename surfaceTensionTable::const_iterator iter =
                stModels_.find(interfacePair(alpha1, alpha2));
            if (iter == stModels_.cend())
            {
                WarningInFunction
                    << "Cannot find interface " << interfacePair(alpha1, alpha2)
                    << " in list of surface tension models" << endl;
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiphaseInterfaceCompressibleSystem::~multiphaseInterfaceCompressibleSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multiphaseInterfaceCompressibleSystem::update()
{
    multiphaseCompressibleSystem::update();
}
//     decode();
//
//     PtrList<surfaceScalarField> rhosOwn(rhos_.size());
//     PtrList<surfaceScalarField> rhosNei(rhos_.size());
//     PtrList<surfaceScalarField> alphasOwn(rhos_.size());
//     PtrList<surfaceScalarField> alphasNei(rhos_.size());
//
//     surfaceScalarField rhoOwn
//     (
//         surfaceScalarField::New
//         (
//             "rhoOwn",
//             mesh(),
//             dimensionedScalar(dimDensity, 0.0)
//         )
//     );
//     surfaceScalarField rhoNei
//     (
//         surfaceScalarField::New
//         (
//             "rhoNei",
//             mesh(),
//             dimensionedScalar(dimDensity, 0.0)
//         )
//     );
//     forAll(rhos_, phasei)
//     {
//         autoPtr<ReconstructionScheme<scalar>> alphaLimiter
//         (
//             ReconstructionScheme<scalar>::New
//             (
//                 alphas_[phasei],
//                 "alpha",
//                 alphas_[phasei].group(),
//                 true
//             )
//         );
//         alphasOwn.set(phasei, alphaLimiter->interpolateOwn());
//         alphasNei.set(phasei, alphaLimiter->interpolateNei());
//
//         autoPtr<ReconstructionScheme<scalar>> rhoLimiter
//         (
//             ReconstructionScheme<scalar>::New
//             (
//                 rhos_[phasei],
//                 "rho",
//                 rhos_[phasei].group(),
//                 true
//             )
//         );
//         rhosOwn.set(phasei, rhoLimiter->interpolateOwn());
//         rhosNei.set(phasei, rhoLimiter->interpolateNei());
//         fluxScheme::correctPhaseFields
//         (
//             alphas_[phasei],
//             rhosOwn[phasei], rhosNei[phasei],
//             thermo_.thermo(phasei).residualAlpha().value()
//         );
//
//         tmp<surfaceScalarField> talphaRhoOwn = surfaceScalarField::New
//         (
//             rhoLimiter->ownName(alphaRhos_[phasei].name()),
//             alphasOwn[phasei]*rhosOwn[phasei]
//         );
//         tmp<surfaceScalarField> talphaRhoNei = surfaceScalarField::New
//         (
//             rhoLimiter->neiName(alphaRhos_[phasei].name()),
//             alphasNei[phasei]*rhosNei[phasei]
//         );
//
//         if (mesh().cacheTemporaryObject(talphaRhoOwn().name()))
//         {
//             mesh().cacheTemporaryObject(talphaRhoOwn.ref());
//             mesh().cacheTemporaryObject(talphaRhoNei.ref());
//         }
//
//         rhoOwn += talphaRhoOwn;
//         rhoNei += talphaRhoNei;
//     }
//     if (mesh().cacheTemporaryObject(rhoOwn.name()))
//     {
//         mesh().cacheTemporaryObject(rhoOwn);
//         mesh().cacheTemporaryObject(rhoNei);
//     }
//
//     fluxScheme_->update
//     (
//         rhoOwn,
//         rhoNei,
//         U_,
//         e_,
//         p_,
//         speedOfSound()(),
//         phi_,
//         rhoPhi_,
//         rhoUPhi_,
//         rhoEPhi_
//     );
//
//
//     // Limit alpha flux
//     volScalarField divPhi(fvc::div(phi_));
//     surfaceScalarField phi(phi_);
//     this->storeAndBlendDelta(phi);
//
//     UPtrList<const volScalarField> alphas(alphas_.size());
//     forAll(alphas_, phasei)
//     {
//         alphas.set(phasei, &alphas_[phasei]);
//         volScalarField alphaOld(alphas_[phasei]);
//         this->blendOld(alphaOld);
//         alphaOld.storeOldTime();
//
//         surfaceScalarField& alphaPhi = alphaPhis_[phasei];
//         alphaPhi =
//             fluxScheme_->flux(alphasOwn[phasei], alphasNei[phasei], phi_);
//         this->storeAndBlendDelta(alphaPhi);
//
//         MULES::limit
//         (
//             1.0/mesh().time().deltaT().value(),
//             geometricOneField(),
//             alphaOld,
//             phi,
//             alphaPhi,
//             zeroField(),
//             zeroField(),//(-divPhi*alphas_[phasei])(),
//             oneField(),
//             zeroField(),
//             false
//         );
//         alphaPhi = this->calcAndStoreDelta(alphaPhi);
//     }
//     MULES::limitSum(alphas, alphaPhis_, phi_);
//
//     // Update phase mass fluxes
//     forAll(alphas_, phasei)
//     {
//         alphaRhoPhis_[phasei] =
//             fluxScheme_->flux
//             (
//                 rhosOwn[phasei],
//                 rhosNei[phasei],
//                 alphaPhis_[phasei]
//             );
//     }
//
//     PtrList<surfaceScalarField> alphaRhoPhiUDs(alphaRhoPhis_.size());
//     forAll(alphaRhoPhis_, phasei)
//     {
//         alphaRhoPhiUDs.set
//         (
//             phasei,
//             upwind<scalar>(mesh(), alphaPhis_[phasei]).flux(rhos_[phasei])
//         );
//
//         alphaRhoPhis_[phasei] -= alphaRhoPhiUDs[phasei];
//     }
//
//     {
//         UPtrList<scalarField> alphaRhoPhisInternal(alphaRhoPhis_.size());
//
//         forAll(alphaRhoPhisInternal, phasei)
//         {
//             alphaRhoPhisInternal.set(phasei, &alphaRhoPhis_[phasei]);
//         }
//
//         MULES::limitSum(alphaRhoPhisInternal);
//     }
//
//     const surfaceScalarField::Boundary& phibf = phi_.boundaryField();
//     forAll(phibf, patchi)
//     {
//         if (phibf[patchi].coupled())
//         {
//             UPtrList<scalarField> alphaRhoPhisPatch(alphaRhoPhis_.size());
//
//             forAll(alphaRhoPhisPatch, phasei)
//             {
//                 alphaRhoPhisPatch.set
//                 (
//                     phasei,
//                     &alphaRhoPhis_[phasei].boundaryFieldRef()[patchi]
//                 );
//             }
//
//             MULES::limitSum(alphaRhoPhisPatch);
//         }
//     }
//
//     forAll(alphaRhoPhis_, phasei)
//     {
//         alphaRhoPhis_[phasei] += alphaRhoPhiUDs[phasei];
//     }
//
//     thermo_.update();
// }


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
