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

#include "UnsTotalLagrangianGeomSolid.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class IncrementalModel>
UnsTotalLagrangianGeomSolid<IncrementalModel>::UnsTotalLagrangianGeomSolid
(
    const word& type,
    fvMesh& mesh,
    const bool isSolid
)
:
    IncrementalModel(type, mesh, nonLinGeom(), isSolid),
    Ff_
    (
        IOobject
        (
            "Ff",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("I", dimless, I)
    ),
    invFf_
    (
        IOobject
        (
            "Ffinv",
            mesh.time().name(),
            mesh
        ),
        inv(Ff_)
    ),
    relFf_
    (
        IOobject
        (
            "relFf",
            mesh.time().name(),
            mesh
        ),
        I + this->gradDDf().T()
    ),
    Jf_
    (
        IOobject
        (
            "Jf",
            mesh.time().name(),
            mesh
        ),
        det(Ff_)
    ),
    relJf_
    (
        IOobject
        (
            "relJf",
            mesh.time().name(),
            mesh
        ),
        det(relFf_)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class IncrementalModel>
void UnsTotalLagrangianGeomSolid<IncrementalModel>::update
(
    const bool correctSigma
)
{
    IncrementalModel::updateDisplacement();

    if (this->incremental())
    {
        // Total deformation gradient
        Ff_ = Ff_.oldTime() + this->gradDDf().T();
    }
    else
    {
        // Total deformation gradient
        Ff_ = I + this->solutionGradDf().T();
    }

    // Jacobian of the deformation gradient
    Jf_ = det(Ff_);

    this->checkEnforceLinear(Jf_);

    if (!this->enforceLinear())
    {
        invFf_ = inv(Ff_);

        // Relative deformation gradient
        relFf_ = Ff_ & invFf_.oldTime();

        // Relative Jacobian (Jacobian of relative deformation gradient)
        relJf_ = det(relFf_);
    }
    else
    {
        Ff_ = Ff_.oldTime();
        Jf_ = Jf_.oldTime();
        invFf_ = invFf_.oldTime();

        relFf_ = tensor::I;
        relJf_ = 1.0;
    }

    if (correctSigma)
    {
        this->mechanical().correct(this->sigmaf());
    }
}


template<class IncrementalModel>
tmp<surfaceTensorField>
UnsTotalLagrangianGeomSolid<IncrementalModel>::Pf() const
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
    const surfaceSymmTensorField& sigmaf = this->sigmaf();

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
        return tPiolaf;
    }

    forAll(Piolaf, facei)
    {
        Piolaf[facei] = Jf_[facei]*(tensor(sigmaf[facei]) & T(invFf_[facei]));
    }
    surfaceTensorField::Boundary& bPiolaf = Piolaf.boundaryFieldRef();
    forAll(bPiolaf, patchi)
    {
        fvsPatchTensorField& pPiolaf = bPiolaf[patchi];
        const fvsPatchSymmTensorField& psigmaf = sigmaf.boundaryField()[patchi];
        const fvsPatchScalarField& pJf = Jf_.boundaryField()[patchi];
        const fvsPatchTensorField& pinvFf = invFf_.boundaryField()[patchi];
        forAll(pPiolaf, facei)
        {
            pPiolaf[facei] =
                pJf[facei]*(tensor(psigmaf[facei]) & T(pinvFf[facei]));
        }
    }

    return tPiolaf;
}


template<class IncrementalModel>
tmp<tensorField>
UnsTotalLagrangianGeomSolid<IncrementalModel>::Pf(const fvPatch& patch) const
{
    const label patchi = patch.index();
    tmp<tensorField> tPiola(new tensorField(patch.size()));
    tensorField& Piola = tPiola.ref();

    const symmTensorField& sigma = this->sigmaf().boundaryField()[patchi];

    if (this->enforceLinear())
    {
        forAll(Piola, facei)
        {
            Piola[facei] = sigma[facei];
        }
        return tPiola;
    }

    const tensorField& invF = this->invFf_.boundaryField()[patchi];
    const scalarField& J = this->Jf_.boundaryField()[patchi];

    forAll(Piola, facei)
    {
        Piola[facei] = J[facei]*(sigma[facei] & T(invF[facei]));
    }
    return tPiola;
}


template<class IncrementalModel>
tmp<vectorField>
UnsTotalLagrangianGeomSolid<IncrementalModel>::nf
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
        tmp<tensorField> pFinvT(inv(Ff_.boundaryField()[patchID])().T());

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
