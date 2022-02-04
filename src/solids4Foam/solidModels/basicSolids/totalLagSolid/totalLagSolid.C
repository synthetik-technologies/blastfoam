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

#include "totalLagSolid.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class IncrementalModel>
void totalLagSolid<IncrementalModel>::update()
{
    IncrementalModel::updateDisplacement();

    if (this->incremental())
    {
        // Total deformation gradient
        relF_ = I + this->gradDD().T();

        // Relative deformation gradient
        F_ = relF_ & F_.oldTime();
    }
    else
    {
        // Total deformation gradient
        F_ = I + this->gradD().T();

        // Relative deformation gradient
        relF_ = F_ & Finv_.oldTime();
    }

    // Inverse of the deformation gradient
    Finv_ = inv(F_);

    // Inverse relative deformation gradient
    relFinv_ = inv(relF_);

    // Jacobian of the deformation gradient
    J_ = det(F_);

    // Relative Jacobian (Jacobian of relative deformation gradient)
    relJ_ = det(relF_);

    // Update stress
    this->mechanical().correct(this->sigma());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class IncrementalModel>
totalLagSolid<IncrementalModel>::totalLagSolid
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
    relF_
    (
        IOobject
        (
            "relF",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        I + this->gradDD().T()
    ),
    Finv_
    (
        IOobject
        (
            "Finv",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        inv(F_)
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
    J_
    (
        IOobject
        (
            "J",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        det(F_)
    ),
    relJ_
    (
        IOobject
        (
            "relJ",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        det(relF_)
    ),
    impK_("impK", this->mechanical().elasticModulus()),
    impKf_("impKf", fvc::interpolate(impK_))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class IncrementalModel>
tmp<vectorField> totalLagSolid<IncrementalModel>::tractionBoundarySnGrad
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
    const tensorField& pGradD =
        IncrementalModel::solutionGradD().boundaryField()[patchID];

    // Patch Cauchy stress
    const symmTensorField& pSigma = this->sigma().boundaryField()[patchID];

    // Patch total deformation gradient inverse
    const tensorField& pFinv = Finv_.boundaryField()[patchID];

    // Patch unit normals (initial configuration)
    const vectorField n(patch.nf());

    // Patch unit normals (deformed configuration)
    vectorField nCurrent(pFinv.T() & n);
    nCurrent /= mag(nCurrent);

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - nCurrent*pressure)
              - (nCurrent & pSigma)
              + pimpK*(n & pGradD)
            )/pimpK
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels
} // End namespace Foam

// ************************************************************************* //
