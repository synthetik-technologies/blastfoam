/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "updatedLagSolid.H"
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
void updatedLagSolid<IncrementalModel>::update
(
    const bool correctSigma
)
{
    if (this->incremental())
    {
        // Update the total displacement
        this->D() = this->D().oldTime() + this->DD();
    }
    else
    {
        // Update the displacement increment
        this->DD() = this->D() - this->D().oldTime();
    }

    // Interpolate DD to pointDD
    this->mechanical().interpolate(this->DD(), this->pointDD(), false);

    // Update gradient of displacement increment
    this->mechanical().grad(this->DD(), this->pointDD(), this->gradDD());

    // Update the gradient of total displacement
    this->gradD() = this->gradD().oldTime() + this->gradDD();

    // Relative deformation gradient
    relF_ = I + this->gradDD().T();

    // Inverse relative deformation gradient
    relFinv_ = inv(relF_);

    // Total deformation gradient
    F_ = relF_ & F_.oldTime();

    // Relative Jacobian
    relJ_ = det(relF_);

    // Jacobian of deformation gradient
    J_ = relJ_*J_.oldTime();

    if (correctSigma)
    {
        this->mechanical().correct(this->sigma());
    }

    impK_ = this->mechanical().impK();
    impKf_ = this->mechanical().impKf();

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class IncrementalModel>
updatedLagSolid<IncrementalModel>::updatedLagSolid
(
    const word& type,
    dynamicFvMesh& mesh,
    const bool isSolid
)
:
    IncrementalModel(type, mesh, nonLinGeom(), isSolid),
    F_
    (
        IOobject
        (
            "F",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("I", dimless, I)
    ),
    J_
    (
        IOobject
        (
            "J",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        det(F_)
    ),
    relF_
    (
        IOobject
        (
            "relF",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        I + this->gradDD().T()
    ),
    relFinv_
    (
        IOobject
        (
            "relFinv",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        inv(relF_)
    ),
    relJ_
    (
        IOobject
        (
            "relJ",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        det(relF_)
    ),
    impK_(this->mechanical().impK()),
    impKf_(this->mechanical().impKf())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class IncrementalModel>
tmp<vectorField> updatedLagSolid<IncrementalModel>::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch implicit stiffness field
    const scalarField& pimpK = impK_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& pgradDD =
        this->solutionGradD().boundaryField()[patchID];

    // Patch stress
    const symmTensorField& psigma =
        this->sigma().boundaryField()[patchID];

    // Patch relative deformation gradient inverse
    const tensorField& prelFinv = relFinv_.boundaryField()[patchID];

    // Patch relative Jacobian
    const scalarField& prelJ = relJ_.boundaryField()[patchID];

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


template<class IncrementalModel>
void updatedLagSolid<IncrementalModel>::updateTotalFields()
{
    // Density
    this->rho() = this->rho().oldTime()/relJ_;

    // Move the mesh to the deformed configuration
    const vectorField oldPoints = this->mesh().points();
    this->moveMesh(oldPoints, this->DD(), this->pointDD());

    IncrementalModel::updateTotalFields();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels
} // End namespace Foam

// ************************************************************************* //
