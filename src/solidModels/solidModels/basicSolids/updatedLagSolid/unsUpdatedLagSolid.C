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

#include "unsUpdatedLagSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "bound.H"
#include "symmetryPolyPatch.H"
#include "twoDPointCorrector.H"
#include "solidTractionFvPatchVectorField.H"
#include "fvcGradf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

template<class IncrementalModel>
void unsUpdatedLagSolid<IncrementalModel>::update
(
    const bool correctSigma
)
{
    updatedLagSolid<IncrementalModel>::update(false);

    // Update gradient of displacement increment
    this->mechanical().grad(this->DD(), this->pointDD(), this->gradDDf_);

    // Update the gradient of total displacement
    this->gradD_ = this->gradD_.oldTime() + this->gradDD_;

    // Relative deformation gradient
    relFf_ = I + this->gradDD_.T();

    // Total deformation gradient
    Ff_ = relFf_ & Ff_.oldTime();

    // Relative Jacobian
    relJf_ = det(relFf_);

    // Jacobian of deformation gradient
    Jf_ = relJf_*Jf_.oldTime();

    if (correctSigma)
    {
        this->mechanical().correct(this->sigmaf_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class IncrementalModel>
unsUpdatedLagSolid<IncrementalModel>::unsUpdatedLagSolid
(
    const word& type,
    dynamicFvMesh& mesh,
    const bool isSolid
)
:
    updatedLagSolid<IncrementalModel>(type, mesh, isSolid),
    Ff_
    (
        IOobject
        (
            "Ff",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(this->F_)
    ),
    relFf_
    (
        IOobject
        (
            "relFf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(this->relF_)
    ),
    Finvf_
    (
        IOobject
        (
            "Finvf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(this->Finv_)
    ),
    relFinvf_
    (
        IOobject
        (
            "relFinvf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(this->relF_)
    ),
    Jf_
    (
        IOobject
        (
            "Jf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(this->J_)
    ),
    relJf_
    (
        IOobject
        (
            "relJf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(this->relJ_)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class IncrementalModel>
tmp<surfaceTensorField> unsUpdatedLagSolid<IncrementalModel>::relPf() const
{
    tmp<surfaceTensorField> tPiolaf
    (
        surfaceTensorField::New
        (
            "Pf",
            this->mesh(),
            dimensionedTensor(this->sigma().dimensions(), Zero)
        )
    );
    surfaceTensorField& Piolaf = tPiolaf.ref();
    const surfaceSymmTensorField& sigmaf = this->sigmaf_;

    if (this->enforceLinear())
    {
        forAll(Piolaf, facei)
        {
            Piolaf[facei] = sigmaf[facei];
        }
        surfaceTensorField::Boundary& bPiolaf = Piolaf.boundaryFieldRef();
        forAll(bPiolaf, patchi)
        {
            fvsPatchTensorField& pPiolaf = bPiolaf[patchi];
            const fvsPatchSymmTensorField& psigmaf = sigmaf.boundaryField()[patchi];
            forAll(pPiolaf, facei)
            {
                pPiolaf[facei] = psigmaf[facei];
            }
        }
        return tPiola;
    }

    forAll(Piolaf, facei)
    {
        Piolaf[facei] = relJf_[facei]*(inv(relFf_[facei]) & tensor(sigmaf[facei]));
    }
    surfaceTensorField::Boundary& bPiolaf = Piolaf.boundaryFieldRef();
    forAll(bPiolaf, patchi)
    {
        fvsPatchTensorField& pPiolaf = bPiolaf[patchi];
        const fvsPatchSymmTensorField& psigmaf = sigmaf.boundaryField()[patchi];
        const fvsPatchScalarField& prelJf = relJf_.boundaryField()[patchi];
        const fvsPatchTensorField& prelFf = relFf_.boundaryField()[patchi];
        forAll(pPiolaf, facei)
        {
            pPiolaf[facei] = prelJf[facei]*(inv(prelFf[facei]) & tensor(psigmaf[facei]));
        }
    }

    return tPiola;
}


template<class IncrementalModel>
tmp<vectorField> unsUpdatedLagSolid<IncrementalModel>::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch implicit stiffness field
    const scalarField& pimpK = this->impKf_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& pgradDD =
        this->solutionGradDf().boundaryField()[patchID];

    // Patch stress
    const symmTensorField& psigma =
        this->sigmaf_.boundaryField()[patchID];


    // Patch unit normals (updated configuration)
    vectorField n(patch.nf());

    if (!this->enforceLinear())
    {
        // Patch relative deformation gradient inverse
        tmp<tensorField> pRelFinvT(inv(relFf_.boundaryField()[patchID])().T());

        // Patch unit normals (deformed configuration)
        n = pRelFinvT & n;
        n /= mag(n);
    }

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - n*pressure)
              - (n & psigma)
              + (patch.nf() & (pimpK*pgradDD))
            )/pimpK
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels
} // End namespace Foam

// ************************************************************************* //
