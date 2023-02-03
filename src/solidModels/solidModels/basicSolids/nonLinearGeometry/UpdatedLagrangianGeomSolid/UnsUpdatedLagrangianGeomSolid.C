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

#include "UnsUpdatedLagrangianGeomSolid.H"
#include "fvcInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class IncrementalModel>
UnsUpdatedLagrangianGeomSolid<IncrementalModel>::UnsUpdatedLagrangianGeomSolid
(
    const word& type,
    dynamicFvMesh& mesh,
    const bool isSolid
)
:
    IncrementalModel(type, mesh, nonLinGeom(), isSolid),
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
        I + this->gradDf().T()
    ),
    relFf_
    (
        IOobject
        (
            "relFf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        I + this->gradDDf().T()
    ),
    invRelFf_
    (
        IOobject
        (
            "invRelFf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        inv(relFf_)
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
        det(Ff_)
    ),
    relJf_
    (
        IOobject
        (
            "relJf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        det(relFf_)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class IncrementalModel>
void UnsUpdatedLagrangianGeomSolid<IncrementalModel>::update
(
    const bool correctSigma
)
{
    IncrementalModel::updateDisplacement();

    // Update gradient of displacement increment
    this->mechanical().grad(this->DD(), this->gradDD());
    this->mechanical().grad(this->DD(), this->pointDD(), this->gradDDf());

    // Update the gradient of total displacement
    this->gradD() = this->gradD().oldTime() + this->gradDD();
    this->gradDf() = this->gradDf().oldTime() + this->gradDDf();

    // Relative deformation gradient
    relFf_ = I + this->gradDDf().T();

    // Total deformation gradient
    Ff_ = relFf_ & Ff_.oldTime();

    // Relative Jacobian
    relJf_ = det(relFf_);

    // Jacobian of deformation gradient
    Jf_ = relJf_*Jf_.oldTime();

    this->checkEnforceLinear(Jf_);

    if (this->enforceLinear())
    {
        // Relative deformation gradient
        relFf_ = tensor::I;

        // Relative Jacobian (Jacobian of relative deformation gradient)
        relJf_ = det(relFf_);

        invRelFf_ = tensor::I;
    }
    else
    {
        invRelFf_ = inv(relFf_);
    }


    if (correctSigma)
    {
        this->mechanical().correct(this->sigmaf());
    }
}


template<class IncrementalModel>
tmp<surfaceTensorField>
UnsUpdatedLagrangianGeomSolid<IncrementalModel>::Pf() const
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
        Piolaf[facei] =
            relJf_[facei]*(tensor(sigmaf[facei]) & T(invRelFf_[facei]));
    }
    surfaceTensorField::Boundary& bPiolaf = Piolaf.boundaryFieldRef();
    forAll(bPiolaf, patchi)
    {
        fvsPatchTensorField& pPiolaf = bPiolaf[patchi];
        const fvsPatchSymmTensorField& psigmaf = sigmaf.boundaryField()[patchi];
        const fvsPatchScalarField& prelJf = relJf_.boundaryField()[patchi];
        const fvsPatchTensorField& pinvRelFf =
            invRelFf_.boundaryField()[patchi];
        forAll(pPiolaf, facei)
        {
            pPiolaf[facei] =
                prelJf[facei]*(tensor(psigmaf[facei]) & T(pinvRelFf[facei]));
        }
    }

    return tPiolaf;
}


template<class IncrementalModel>
tmp<tensorField>
UnsUpdatedLagrangianGeomSolid<IncrementalModel>::Pf(const fvPatch& patch) const
{
    const label patchi = patch.index();
    tmp<tensorField> tpPiola(new tensorField(patch.size()));
    tensorField& pPiola = tpPiola.ref();
    const symmTensorField& psigma = this->sigmaf().boundaryField()[patchi];
    const tensorField& prelF = this->relFf_.boundaryField()[patchi];
    const scalarField& prelJ = this->relJf_.boundaryField()[patchi];

    if (this->enforceLinear())
    {
        forAll(pPiola, facei)
        {
            pPiola[facei] = psigma[facei];
        }
        return tpPiola;
    }

    forAll(pPiola, facei)
    {
        pPiola[facei] = prelJ[facei]*(psigma[facei] & T(inv(prelF[facei])));
    }
    return tpPiola;
}



template<class IncrementalModel>
tmp<vectorField>
UnsUpdatedLagrangianGeomSolid<IncrementalModel>::nf(const fvPatch& patch) const
{
    tmp<vectorField> n(patch.nf());
    if (!this->enforceLinear())
    {
        // Patch index
        const label patchID = patch.index();

        // Patch relative deformation gradient inverse
        tmp<tensorField> pRelFinvT(invRelFf_.boundaryField()[patchID].T());

        // Patch unit normals (deformed configuration)
        n.ref() = pRelFinvT & n();
        n.ref() /= mag(n());
    }
    return n;
}


template<class IncrementalModel>
void UnsUpdatedLagrangianGeomSolid<IncrementalModel>::updateTotalFields()
{
    // Density
    this->rho() = this->rho().oldTime()/fvc::surfVolInterpolate(relJf_);

    // Move the mesh to the deformed configuration
    const vectorField oldPoints = this->mesh().points();
    this->moveMesh(oldPoints, this->DD(), this->pointDD());

    IncrementalModel::updateTotalFields();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels
} // End namespace Foam

// ************************************************************************* //
