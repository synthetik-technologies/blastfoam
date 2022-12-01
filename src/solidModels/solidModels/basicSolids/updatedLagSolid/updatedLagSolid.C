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
#include "globalPolyBoundaryMesh.H"

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
    IncrementalModel::updateDisplacement();

    // Relative deformation gradient
    relF_ = I + this->gradDD().T();

    // Total deformation gradient
    F_ = relF_ & F_.oldTime();

    // Relative Jacobian
    relJ_ = det(relF_);

    // Jacobian of deformation gradient
    J_ = relJ_*J_.oldTime();

    this->checkEnforceLinear(J_);

    if (correctSigma)
    {
        this->mechanical().correct(this->sigma());
        impK_ = this->mechanical().impK();
        impKf_ = this->mechanical().impKf();
    }
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
            IOobject::AUTO_WRITE
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
{
    this->globalPatches().setDisplacementField(this->mesh().name(), "pointD");
    this->globalPatches().setInverseDisplacement(this->mesh().name(), true);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class IncrementalModel>
tmp<volTensorField> updatedLagSolid<IncrementalModel>::relP() const
{
    tmp<volTensorField> tPiola
    (
        volTensorField::New
        (
            "P",
            this->mesh(),
            dimensionedTensor(this->sigma().dimensions(), Zero)
        )
    );
    volTensorField& Piola = tPiola.ref();
    const volSymmTensorField& sigma = this->sigma();

    if (this->enforceLinear())
    {
        forAll(Piola, celli)
        {
            Piola[celli] = sigma[celli];
        }
        volTensorField::Boundary& bPiola = Piola.boundaryFieldRef();
        forAll(bPiola, patchi)
        {
            fvPatchTensorField& pPiola = bPiola[patchi];
            const fvPatchSymmTensorField& psigma = sigma.boundaryField()[patchi];
            forAll(pPiola, facei)
            {
                pPiola[facei] = psigma[facei];
            }
        }
        return tPiola;
    }

    forAll(Piola, celli)
    {
        Piola[celli] = relJ_[celli]*(inv(relF_[celli]) & tensor(sigma[celli]));
    }
    volTensorField::Boundary& bPiola = Piola.boundaryFieldRef();
    forAll(bPiola, patchi)
    {
        fvPatchTensorField& pPiola = bPiola[patchi];
        const fvPatchSymmTensorField& psigma = sigma.boundaryField()[patchi];
        const fvPatchScalarField& prelJ = relJ_.boundaryField()[patchi];
        const fvPatchTensorField& prelF = relF_.boundaryField()[patchi];
        forAll(pPiola, facei)
        {
            pPiola[facei] = prelJ[facei]*(inv(prelF[facei]) & tensor(psigma[facei]));
        }
    }

    return tPiola;
}


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

    // Patch unit normals (updated configuration)
    vectorField n(patch.nf());

    if (!this->enforceLinear())
    {
        // Patch relative deformation gradient inverse
        tmp<tensorField> pRelFinvT(inv(relF_.boundaryField()[patchID])().T());

        // Patch unit normals (deformed configuration)
        n = pRelFinvT & n;
    }

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - n*pressure)
              - (n & psigma)
            )/pimpK
          + (patch.nf() & pgradDD)
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
