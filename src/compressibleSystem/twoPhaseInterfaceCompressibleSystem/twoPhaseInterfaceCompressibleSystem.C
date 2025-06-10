/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2025
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

#include "twoPhaseInterfaceCompressibleSystem.H"
#include "addToRunTimeSelectionTable.H"
#include "EulerDdtScheme.H"
#include "gaussConvectionScheme.H"
#include "contactAngleFvPatchScalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(twoPhaseInterfaceCompressibleSystem, 0);
    addToRunTimeSelectionTable
    (
        compressibleSystem,
        twoPhaseInterfaceCompressibleSystem,
        twoPhase
    );
}


// * * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * //

void Foam::twoPhaseInterfaceCompressibleSystem::addSources
(
    volVectorField::Internal& rhoUSource,
    volScalarField::Internal& rhoESource
) const
{

    twoPhaseCompressibleSystem::addSources(rhoUSource, rhoESource);

    tmp<volVectorField> stF
    (
        -surfaceTensionModel_->sigma()*fvc::div(nHatf_)*fvc::grad(alpha1_)
    );
    // (
    //     fvc::reconstruct
    //     (
    //       - fvc::interpolate
    //         (
    //             surfaceTensionModel_->sigma()
    //            *fvc::div(nHatf_)
    //         )
    //        *fvc::snGrad(alpha1_)
    //        *mesh().magSf()
    //     )
    // );
    rhoUSource -= stF();
    rhoESource -= stF() & U_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseInterfaceCompressibleSystem::twoPhaseInterfaceCompressibleSystem
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    twoPhaseCompressibleSystem(dict, mesh),
    surfaceTensionModel_(surfaceTensionModel::New(dict, mesh)),

    deltaN_
    (
        "deltaN",
        1e-8/cbrt(min(this->mesh().V()))
    ),
    nHatf_
    (
        IOobject
        (
            "nHatf",
            mesh.time().name(),
            this->mesh()
        ),
        mesh,
        dimensionedScalar(dimArea, 0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseInterfaceCompressibleSystem::~twoPhaseInterfaceCompressibleSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::twoPhaseInterfaceCompressibleSystem::update()
{
    // Calculate fluxes
    twoPhaseCompressibleSystem::update();

    // Calculate interface normal
    const volVectorField gradAlpha(fvc::grad(alpha1_));
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));

    // Correct contact angle
    {
        const fvBoundaryMesh& bmesh = this->mesh().boundary();
        const volVectorField::Boundary& bU = U_.boundaryField();
        const surfaceVectorField::Boundary& bgradAlpha =
            gradAlphaf.boundaryField();

        surfaceVectorField::Boundary& bnHatfv = nHatfv.boundaryFieldRef();
        volScalarField::Boundary& balpha1 = alpha1_.boundaryFieldRef();
        volScalarField::Boundary& balpha2 = alpha2_.boundaryFieldRef();
        forAll(balpha1, patchi)
        {
            if (isA<contactAngleFvPatchScalarField>(balpha1[patchi]))
            {
                contactAngleFvPatchScalarField& palpha1 =
                    refCast<contactAngleFvPatchScalarField>
                    (
                        balpha1[patchi]
                    );
                fvsPatchVectorField& pnHat = bnHatfv[patchi];
                const scalarField cosTheta
                (
                    palpha1.cosTheta(bU[patchi], pnHat)
                );
                const vectorField nf(bmesh[patchi].nf());
                scalarField a(pnHat.size());
                scalarField b(pnHat.size());
                forAll(a, facei)
                {
                    const scalar a12 = pnHat[facei] & nf[facei];
                    const scalar b1 = cosTheta[facei];
                    const scalar b2 = cos(acos(a12) - acos(cosTheta[facei]));
                    const scalar det = 1.0 - sqr(a12);

                    const scalar a = (b1 - a12*b2)/det;
                    const scalar b = (b2 - a12*b1)/det;

                    pnHat[facei] = a*nf[facei] + b*pnHat[facei];
                    pnHat[facei] /= mag(pnHat[facei]) + deltaN_.value();
                }
                palpha1.gradient() = (nf & pnHat)*mag(bgradAlpha[patchi]);
                palpha1.evaluate();
                balpha2[patchi] = palpha1;
            }
        }
    }

    nHatf_ = nHatfv & mesh().Sf();
}

// ************************************************************************* //
