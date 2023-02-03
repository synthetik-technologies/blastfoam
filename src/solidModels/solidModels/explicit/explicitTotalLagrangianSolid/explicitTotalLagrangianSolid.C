/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "explicitTotalLagrangianSolid.H"
#include "ReconstructionScheme.H"
#include "wedgePolyPatch.H"
#include "symmetryPolyPatch.H"
#include "symmetryPlanePolyPatch.H"
#include "solidTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(explicitTotalLagrangianSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, explicitTotalLagrangianSolid, dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

explicitTotalLagrangianSolid::explicitTotalLagrangianSolid
(
    dynamicFvMesh& mesh
)
:
    TotalLagrangianGeomSolid<explicitNonLinearSolid>(typeName, mesh)
{
    // Update stress
    update();

    // Update initial acceleration
    a_ = fvc::div(sigma(), "div(sigma)")/rho();
    a_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void explicitTotalLagrangianSolid::update(const bool correctSigma)
{
    TotalLagrangianGeomSolid<explicitNonLinearSolid>::update(correctSigma);

    if (correctSigma)
    {
        this->impKf_ = this->mechanical().impKf();
        wavespeed_ = sqrt(this->impKf_/fvc::interpolate(rho()));
    }
    this->updateWavespeeds();
}


tmp<volVectorField> explicitTotalLagrangianSolid::divStress() const
{
    return fvc::div(this->P(), "div(sigma)");
}


tmp<surfaceVectorField> explicitTotalLagrangianSolid::tractionSf() const
{
    // Surface normals
    // const surfaceVectorField N(this->mesh().Sf()/this->mesh().magSf());
    // surfaceVectorField n(fvc::interpolate(T(inv(this->F_))) & N);
    // n /= mag(n);
    // surfaceTensorField nn(n*n);
    // surfaceTensorField stabRhoU
    // (
    //     wavespeed_*nn
    //   + sWavespeed_*(I - nn)
    // );
    //
    // volTensorField P(this->P());
    //
    // // Reconstruction of Piola tensor
    // autoPtr<ReconstructionScheme<tensor>> PLimiter
    // (
    //     ReconstructionScheme<tensor>::New(P, "P")
    // );
    // surfaceVectorField tractionOwn(PLimiter->interpolateOwn() & N);
    // surfaceVectorField tractionNei(PLimiter->interpolateNei() & N);
    //
    // // Reconstruction of Piola tensor
    // autoPtr<ReconstructionScheme<vector>> ULimiter
    // (
    //     ReconstructionScheme<vector>::New(this->U(), "U")
    // );
    // surfaceVectorField UOwn(ULimiter->interpolateOwn());
    // surfaceVectorField UNei(ULimiter->interpolateNei());
    //
    //
    // // Acoustic Riemann solver
    // surfaceVectorField tractionC
    // (
    //     0.5
    //    *(
    //        tractionOwn + tractionNei
    //      + fvc::interpolate(this->rho())*(stabRhoU & (UNei - UOwn))
    //     )
    // );
    //
    // const volScalarField::Boundary& brho(this->rho().boundaryField());
    // const volVectorField::Boundary& bU(this->U().boundaryField());
    // surfaceVectorField::Boundary& btractionC(tractionC.boundaryFieldRef());
    //
    //
    // forAll(btractionC, patchi)
    // {
    //     const polyPatch& p = mesh().boundaryMesh()[patchi];
    //     const fvPatchField<vector>& pDD(this->DD().boundaryField()[patchi]);
    //
    //     if (isA<tractionBase>(pDD))
    //     {
    //         const vectorField pn(this->nf(mesh().boundary()[patchi]));
    //         const tractionBase& tb = dynamicCast<const tractionBase>(pDD);
    //         vectorField tp((tb.traction() - pn*tb.pressure()));
    //
    //         btractionC[patchi] = tp;
    //     }
    //     else if (pDD.fixesValue())
    //     {
    //         btractionC[patchi] ==
    //             tractionOwn.boundaryField()[patchi]
    //           + brho[patchi]
    //            *(
    //                 stabRhoU.boundaryField()[patchi]
    //               & (
    //                     pDD/mesh().time().deltaTValue()
    //                   - UOwn.boundaryField()[patchi]
    //                 )
    //             );
    //     }
    //     else if
    //     (
    //         isA<symmetryPolyPatch>(p)
    //      || isA<symmetryPlanePolyPatch>(p)
    //     )
    //     {
    //         const vectorField pn(this->nf(mesh().boundary()[patchi]));
    //         btractionC[patchi] =
    //             (pn*pn)
    //           & (
    //                 tractionOwn.boundaryField()[patchi]
    //               - wavespeed_.boundaryField()[patchi]
    //                *brho[patchi]*UOwn.boundaryField()[patchi]
    //             );
    //     }
    //     else if (!pDD.coupled())
    //     {
    //         btractionC[patchi] = tractionOwn.boundaryField()[patchi];
    //     }
    // }
    // return tractionC*mesh().magSf();

    return fvc::dotInterpolate(this->mesh().Sf(), this->P());
}


bool explicitTotalLagrangianSolid::evolve()
{
    this->solveMomentum();
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels
} // End namespace Foam

// ************************************************************************* //
