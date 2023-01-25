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

#include "TotalLagrangianGeomSolid.H"
#include "transformGeometricField.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class IncrementalModel>
TotalLagrangianGeomSolid<IncrementalModel>::TotalLagrangianGeomSolid
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
            mesh
        ),
        I + this->gradDD().T()
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
            mesh
        ),
        det(relF_)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class IncrementalModel>
void TotalLagrangianGeomSolid<IncrementalModel>::update(const bool correctSigma)
{
    IncrementalModel::updateDisplacement();

    if (this->incremental())
    {
        // Total deformation gradient
        F_ = F_.oldTime() + this->gradDD().T();
    }
    else
    {
        // Total deformation gradient
        F_ = I + this->gradD().T();
    }
    F_.correctBoundaryConditions();

    // Jacobian of the deformation gradient
    J_ = det(F_);

    this->checkEnforceLinear(J_);

    if (!this->enforceLinear())
    {
        // Relative deformation gradient
        relF_ = F_ & inv(F_.oldTime());
        relF_.correctBoundaryConditions();

        // Relative Jacobian (Jacobian of relative deformation gradient)
        relJ_ = det(relF_);
    }


    // Update stress
    if (correctSigma)
    {
        this->mechanical().correct(this->sigma());
        this->sigma().correctBoundaryConditions();
    }
}


template<class IncrementalModel>
tmp<volTensorField> TotalLagrangianGeomSolid<IncrementalModel>::P() const
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
        Piola[celli] = J_[celli]*(sigma[celli] & T(inv(F_[celli])));
    }
    volTensorField::Boundary& bPiola = Piola.boundaryFieldRef();
    forAll(bPiola, patchi)
    {
        fvPatchTensorField& pPiola = bPiola[patchi];
        const fvPatchSymmTensorField& psigma = sigma.boundaryField()[patchi];
        const fvPatchScalarField& pJ = J_.boundaryField()[patchi];
        const fvPatchTensorField& pF = F_.boundaryField()[patchi];
        forAll(pPiola, facei)
        {
            pPiola[facei] = pJ[facei]*(psigma[facei] & T(inv(pF[facei])));
        }
    }

    return tPiola;
}


template<class IncrementalModel>
tmp<tensorField>
TotalLagrangianGeomSolid<IncrementalModel>::P(const fvPatch& patch) const
{
    const label patchi = patch.index();
    tmp<tensorField> tPiola(new tensorField(patch.size()));
    tensorField& Piola = tPiola.ref();
    const symmTensorField& sigma = this->sigma().boundaryField()[patchi];
    const tensorField& F = this->F_.boundaryField()[patchi];
    const scalarField& J = this->J_.boundaryField()[patchi];

    if (this->enforceLinear())
    {
        forAll(Piola, facei)
        {
            Piola[facei] = sigma[facei];
        }
        return tPiola;
    }

    forAll(Piola, facei)
    {
        Piola[facei] = J[facei]*(sigma[facei] & T(inv(F[facei])));
    }
    return tPiola;
}


template<class IncrementalModel>
tmp<vectorField> TotalLagrangianGeomSolid<IncrementalModel>::nf
(
    const fvPatch& patch
) const
{
    tmp<vectorField> n(patch.nf());
    if (!this->enforceLinear())
    {
        // Patch index
        const label patchID = patch.index();

        // Patch relative deformation gradient inverse
        tmp<tensorField> pFinvT(inv(F_.boundaryField()[patchID])().T());

        // Patch unit normals (deformed configuration)
        n.ref() = pFinvT & n();
        n.ref() /= mag(n());
    }
    return n;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels
} // End namespace Foam

// ************************************************************************* //
