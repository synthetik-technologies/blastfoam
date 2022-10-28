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

#include "unsTotalLagSolid.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class IncrementalModel>
void unsTotalLagSolid<IncrementalModel>::update
(
    const bool correctSigma
)
{
    totalLagSolid<IncrementalModel>::update();

    if (this->incremental())
    {
        // Total deformation gradient
        Ff_ = Ff_.oldTime() + this->solutionGradDf().T();
    }
    else
    {
        // Total deformation gradient
        Ff_ = I + this->solutionGradDf().T();
    }

    // Jacobian of the deformation gradient
    Jf_ = det(Ff_);

    // Relative deformation gradient
    relFf_ = Ff_ & inv(Ff_.oldTime());

    // Relative Jacobian (Jacobian of relative deformation gradient)
    relJf_ = det(relFf_);

    if (correctSigma)
    {
        this->mechanical().correct(this->sigmaf_);
        this->impK_ = this->mechanical().impK();
        this->impKf_ = this->mechanical().impKf();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class IncrementalModel>
unsTotalLagSolid<IncrementalModel>::unsTotalLagSolid
(
    const word& type,
    dynamicFvMesh& mesh,
    const bool isSolid
)
:
    totalLagSolid<IncrementalModel>(type, mesh, isSolid),
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
tmp<surfaceTensorField> unsTotalLagSolid<IncrementalModel>::Pf() const
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
        return tPiolaf;
    }

    forAll(Piolaf, facei)
    {
        Piolaf[facei] = Jf_[facei]*(inv(Ff_[facei]) & tensor(sigmaf[facei]));
    }
    surfaceTensorField::Boundary& bPiolaf = Piolaf.boundaryFieldRef();
    forAll(bPiolaf, patchi)
    {
        fvsPatchTensorField& pPiolaf = bPiolaf[patchi];
        const fvsPatchSymmTensorField& psigmaf = sigmaf.boundaryField()[patchi];
        const fvsPatchScalarField& pJf = Jf_.boundaryField()[patchi];
        const fvsPatchTensorField& pFf = Ff_.boundaryField()[patchi];
        forAll(pPiolaf, facei)
        {
            pPiolaf[facei] = pJf[facei]*(inv(pFf[facei]) & tensor(psigmaf[facei]));
        }
    }

    return tPiolaf;
}


template<class IncrementalModel>
tmp<vectorField> unsTotalLagSolid<IncrementalModel>::tractionBoundarySnGrad
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
    const tensorField& pGradD =
        totalLagSolid<IncrementalModel>::solutionGradDf().boundaryField()[patchID];

    // Patch Cauchy stress
    const symmTensorField& pSigma = this->sigmaf_.boundaryField()[patchID];

    // Patch unit normals (initial configuration)
    vectorField n(patch.nf());

    if (!this->enforceLinear())
    {
        // Patch total deformation gradient inverse
        tmp<tensorField> pFinvT(inv(Ff_.boundaryField()[patchID])().T());

        // Patch unit normals (deformed configuration)
        n = pFinvT & n;
        n /= mag(n);
    }

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - n*pressure)
              - (n & pSigma)
              + pimpK*(n & pGradD)
            )/pimpK
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels
} // End namespace Foam

// ************************************************************************* //
