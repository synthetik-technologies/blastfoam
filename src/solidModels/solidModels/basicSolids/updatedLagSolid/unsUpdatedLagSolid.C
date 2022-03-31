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
void unsUpdatedLagSolid<IncrementalModel>::update()
{
    updatedLagSolid<IncrementalModel>::update();

    // Update gradient of displacement increment
    this->mechanical().grad(this->DD(), this->pointDD(), this->gradDDf_);

    // Update the gradient of total displacement
    this->gradD_ = this->gradD_.oldTime() + this->gradDD_;

    // Relative deformation gradient
    relFf_ = I + this->gradDD_.T();

    // Inverse relative deformation gradient
    relFinvf_ = inv(relFf_);

    // Total deformation gradient
    Ff_ = relFf_ & Ff_.oldTime();

    // Relative Jacobian
    relJf_ = det(relFf_);

    // Jacobian of deformation gradient
    Jf_ = relJf_*Jf_.oldTime();

    this->mechanical().correct(this->sigmaf_);
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

    // Patch relative deformation gradient inverse
    const tensorField& prelFinv = relFinvf_.boundaryField()[patchID];

    // Patch relative Jacobian
    const scalarField& prelJ = relJf_.boundaryField()[patchID];

    // Patch unit normals (updated configuration)
    const vectorField n(patch.nf());

    // Patch unit normals (deformed configuration)
    const vectorField nCurrent(prelJ*prelFinv.T() & n);

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - n*pressure)
              - (nCurrent & psigma)
              + (n & (pimpK*pgradDD))
            )/pimpK
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels
} // End namespace Foam

// ************************************************************************* //
