/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "gradientSchemes.H"
#include "fvc.H"
#include "leastSquaresGrad.H"
#include "ReconstructionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(gradientSchemes, 0);


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

gradientSchemes::gradientSchemes
(
    const fvMesh& mesh
)
:
    DemandDrivenMeshObject
    <
        fvMesh,
        MoveableMeshObject,
        gradientSchemes
    >(mesh),
    ops_(mesh),
    AinvPtr_(nullptr),
    AinvLocalPtr_(nullptr)
{}


bool gradientSchemes::movePoints()
{
    deleteDemandDrivenData(AinvPtr_);
    deleteDemandDrivenData(AinvLocalPtr_);
    return true;
}



// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

gradientSchemes::~gradientSchemes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const tensorField& gradientSchemes::distanceMatrix() const
{
    if (!AinvPtr_)
    {
        const fvMesh& mesh = this->mesh();

        AinvPtr_ = new tensorField(mesh.nCells(), Zero);
        tensorField& Ainv = *AinvPtr_;

        const volVectorField& C = mesh.C();
        const labelList& owner = mesh.owner();
        const labelList& neighbour = mesh.neighbour();
        forAll(owner, facei)
        {
            const label own = owner[facei];
            const label nei = neighbour[facei];
            const vector dOwn = C[nei] - C[own];
            const vector dNei  = C[own] - C[nei];

            Ainv[own] += dOwn*dOwn;
            Ainv[nei] += dNei*dNei;
        }

        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& patch = mesh.boundary()[patchi];
            const vectorField pd(patch.delta());

            if (patch.coupled())
            {
                forAll(mesh.boundary()[patchi], facei)
                {
                    const label celli =
                        mesh.boundaryMesh()[patchi].faceCells()[facei];
                    Ainv[celli] += pd[facei]*pd[facei];
                }
            }
        }

        forAll(Ainv, celli)
        {
            Ainv[celli] = stabInv(Ainv[celli]);
        }
    }
    return *AinvPtr_;
}



void gradientSchemes::calcDistanceMatrixLocal() const
{
    if (AinvLocalPtr_ || AinvLocalFixedPtr_)
    {
        FatalErrorInFunction
            << "AinvLocal already calculated" << endl
            << abort(FatalError);
    }

    const fvMesh& mesh = this->mesh();

    AinvLocalPtr_ = new tensorField(mesh.nCells(), Zero);
    AinvLocalFixedPtr_ = new tensorField(mesh.nCells(), Zero);

    tensorField& AinvLocal = *AinvLocalPtr_;
    tensorField& AinvLocalFixed = *AinvLocalFixedPtr_;

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();
    const pointField& points = mesh.points();
    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();
    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];
        const vector dOwn = Cf[facei] - C[own];
        const vector dNei = Cf[facei] - C[nei];

        AinvLocal[own] += dOwn*dOwn;
        AinvLocal[nei] += dNei*dNei;
    }

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        const vectorField pd(patch.fvPatch::delta());

        forAll(mesh.boundary()[patchi], facei)
        {
            const label celli =
                mesh.boundaryMesh()[patchi].faceCells()[facei];

            AinvLocal[celli] += pd[facei]*pd[facei];

            // Fixed value correction
            const label mfacei =
                mesh.boundary()[patchi].start() + facei;
            forAll(mesh.faces()[mfacei], nodei)
            {
                const label pti = mesh.faces()[mfacei][nodei];

                vector d = points[pti] - C[celli];
                AinvLocalFixed[celli] += d*d;

                for (label i = 0; i < 7; i++)
                {
                    scalar si(i);
                    d =
                        (
                            (
                                (si + 1.0)*points[pti]
                              + (
                                    (7.0 - si)
                                   *Cf.boundaryField()[patchi][facei]
                                )
                            )
                        )/8.0 - C[celli];
                    AinvLocalFixed[celli] += d*d;
                }
            }
        }
    }

    AinvLocalFixed += AinvLocal;
    forAll(AinvLocal, celli)
    {
        AinvLocal[celli] = stabInv(AinvLocal[celli]);
        AinvLocalFixed[celli] = stabInv(AinvLocalFixed[celli]);
    }
}


const tensorField& gradientSchemes::distanceMatrixLocal() const
{
    if (!AinvLocalPtr_)
    {
       calcDistanceMatrixLocal();
    }
    return *AinvLocalPtr_;
}


const tensorField& gradientSchemes::distanceMatrixLocalFixed() const
{
    if (!AinvLocalFixedPtr_)
    {
       calcDistanceMatrixLocal();
    }
    return *AinvLocalFixedPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volVectorField> gradientSchemes::gradient
(
    const GeometricField<scalar, fvPatchField, volMesh>& U
)   const
{
    const fvMesh& mesh = this->mesh();

    tmp<volVectorField> tgradU
    (
        volVectorField::New
        (
            "grad(" + U.name() + ")",
            mesh,
            dimensioned<vector>
            (
                "0",
                U.dimensions()/dimLength,
                pTraits<vector>::zero
            )
        )
    );
    volVectorField& gradU = tgradU.ref();

    const volVectorField& C = mesh.C();
    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();
    const tensorField& Ainv = distanceMatrix();
    forAll(owner, faceID)
    {
        const label own = owner[faceID];
        const label nei = neighbour[faceID];
        const vector dOwn = C[nei] - C[own];
        const vector dNei = C[own] - C[nei];

        gradU[own] += Ainv[own] & ((U[nei] - U[own])*dOwn);
        gradU[nei] += Ainv[nei] & ((U[own] - U[nei])*dNei);
    }

    const volScalarField::Boundary& pU(U.boundaryField());
    volVectorField::Boundary& pgradU = gradU.boundaryFieldRef();
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];

        if (mesh.boundary()[patchi].coupled())
        {
            const vectorField pd(mesh.boundary()[patchi].delta());
            const scalarField UNei
            (
                pU[patchi].patchNeighbourField()
            );

            pgradU[patchi] =
                patch.deltaCoeffs()
               *(UNei - pU[patchi].patchInternalField())*patch.nf();

            forAll(mesh.boundary()[patchi], facei)
            {
                const label& celli =
                    mesh.boundaryMesh()[patchi].faceCells()[facei];

                gradU[celli] +=
                    Ainv[celli] & ((UNei[facei] - U[celli])*pd[facei]);
            }
        }
        else
        {
            pgradU[patchi] =
                patch.deltaCoeffs()
               *(pU[patchi] - pU[patchi].patchInternalField())
               *patch.nf();
        }
    }

    return tgradU;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volTensorField> gradientSchemes::gradient
(
    const volVectorField& U
)   const
{
    return fvc::grad(U);
    volVectorField gradUx(gradientSchemes::gradient(U.component(0)));
    volVectorField gradUy(gradientSchemes::gradient(U.component(1)));
    volVectorField gradUz(gradientSchemes::gradient(U.component(2)));

    tmp<volTensorField > tgradU
    (
        volTensorField::New
        (
            "grad(" + U.name() + ")",
            mesh(),
            dimensioned<tensor>
            (
                "0",
                U.dimensions()/dimLength,
                pTraits<tensor>::zero
            )
        )
    );
    volTensorField& gradU = tgradU.ref();

    forAll(gradU, celli)
    {
        gradU[celli] = tensor(gradUx[celli], gradUy[celli], gradUz[celli]);
    }

    volTensorField::Boundary& pgradU(gradU.boundaryFieldRef());
    forAll(pgradU, patchi)
    {
        forAll(pgradU[patchi], facei)
        {
            pgradU[patchi][facei] =
                tensor
                (
                    gradUx.boundaryField()[patchi][facei],
                    gradUy.boundaryField()[patchi][facei],
                    gradUz.boundaryField()[patchi][facei]
                );
        }
    }

    return tgradU;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::gradient
(
    const volTensorField& U,
    volTensorField& UgradX,
    volTensorField& UgradY,
    volTensorField& UgradZ
)   const
{
    UgradX = gradientSchemes::gradient(ops_.decomposeTensorX(U));
    UgradY = gradientSchemes::gradient(ops_.decomposeTensorY(U));
    UgradZ = gradientSchemes::gradient(ops_.decomposeTensorZ(U));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volTensorField> gradientSchemes::localGradient
(
    const volVectorField& U,
    const GeometricField<vector, fvsPatchField, surfaceMesh>& Uf,
    const pointVectorField& pointU
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<volTensorField > tgradU
    (
        volTensorField::New
        (
            "grad(" + U.name() + ")",
            mesh,
            dimensioned<tensor>
            (
                "0",
                U.dimensions()/dimLength,
                Zero
            )
        )
    );
    volTensorField& gradU = tgradU.ref();

    volVectorField gradUx
    (
        volVectorField::New
        (
            "grad" + U.name() + "x",
            mesh,
            dimensionedVector(gradU.dimensions(), Zero)
        )
    );
    volVectorField gradUy
    (
        volVectorField::New
        (
            "grad" + U.name() + "y",
            mesh,
            dimensionedVector(gradU.dimensions(), Zero)
        )
    );
    volVectorField gradUz
    (
        volVectorField::New
        (
            "grad" + U.name() + "z",
            mesh,
            dimensionedVector(gradU.dimensions(), Zero)
        )
    );

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();
    const pointField& points = mesh.points();
    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();
    const tensorField& AinvLocal = distanceMatrixLocal();
    const tensorField& AinvLocalFixed = distanceMatrixLocalFixed();
    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];
        const vector dOwn = Cf[facei] - C[own];
        const vector dNei = Cf[facei] - C[nei];

        gradUx[own] +=
            AinvLocal[own] & ((Uf[facei].x() - U[own].x())*dOwn);
        gradUy[own] +=
            AinvLocal[own] & ((Uf[facei].y() - U[own].y())*dOwn);
        gradUz[own] +=
            AinvLocal[own] & ((Uf[facei].z() - U[own].z())*dOwn);

        gradUx[nei] +=
            AinvLocal[nei] & ((Uf[facei].x() - U[nei].x())*dNei);
        gradUy[nei] +=
            AinvLocal[nei] & ((Uf[facei].y() - U[nei].y())*dNei);
        gradUz[nei] +=
            AinvLocal[nei] & ((Uf[facei].z() - U[nei].z())*dNei);
    }

    // const volVectorField::Boundary& pU(U.boundaryField());
    const surfaceVectorField::Boundary& pUf(Uf.boundaryField());
    // volVectorField::Boundary& pgradUx(gradUx.boundaryFieldRef());
    // volVectorField::Boundary& pgradUy(gradUy.boundaryFieldRef());
    // volVectorField::Boundary& pgradUz(gradUz.boundaryFieldRef());

    forAll(U.boundaryField(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        const polyPatch& ppatch = mesh.boundaryMesh()[patchi];
        // tensorField pgradU
        // (
        //     patch.deltaCoeffs()
        //    *patch.nf()
        //    *(pUf[patchi] - pU[patchi].patchInternalField())
        // );
        // ops_.decomposeTensor
        // (
        //     pgradU,
        //     pgradUx[patchi],
        //     pgradUy[patchi],
        //     pgradUz[patchi]
        // );

        if (pointU.boundaryField()[patchi].fixesValue())
        {
            forAll(C.boundaryField()[patchi], facei)
            {
                const label celli = patch.faceCells()[facei];
                vector d = Cf.boundaryField()[patchi][facei] - C[celli];

                gradUx[celli] +=
                    AinvLocal[celli]
                  & ((pUf[patchi][facei].x() - U[celli].x())*d);

                gradUy[celli] +=
                    AinvLocal[celli]
                  & ((pUf[patchi][facei].y() - U[celli].y())*d);

                gradUz[celli] +=
                    AinvLocal[celli]
                  & ((pUf[patchi][facei].z() - U[celli].z())*d);

                const face& f = ppatch[facei];
                forAll(f, pti)
                {
                    const label pointi = f[pti];
                    vector d = points[pointi] - C[celli];

                    gradUx[celli] +=
                        AinvLocal[celli]
                      & ((pointU[pointi].x() - U[celli].x())*d);

                    gradUy[celli] +=
                        AinvLocal[celli]
                      & ((pointU[pointi].y() - U[celli].y())*d);

                    gradUz[celli] +=
                        AinvLocal[celli]
                      & ((pointU[pointi].z() - U[celli].z())*d);

                    for (label i = 0; i < 7; i++)
                    {
                        scalar si(i);
                        d =
                            (
                                (
                                    (si + 1.0)*points[pointi]
                                  + (
                                        (7.0 - si)
                                       *Cf.boundaryField()[patchi][facei]
                                    )
                                )/8.0
                            ) - C[celli];

                        gradUx[celli] +=
                            AinvLocalFixed[celli]
                          & ((pointU[pointi].x() - U[celli].x())*d);

                        gradUy[celli] +=
                            AinvLocalFixed[celli]
                          & ((pointU[pointi].y() - U[celli].y())*d);

                        gradUz[celli] +=
                            AinvLocalFixed[celli]
                          & ((pointU[pointi].z() - U[celli].z())*d);
                    }
                }
            }
        }
        else
        {
            forAll(C.boundaryField()[patchi], facei)
            {
                const label celli = patch.faceCells()[facei];
                vector d = Cf.boundaryField()[patchi][facei] - C[celli];

                gradUx[celli] +=
                    AinvLocal[celli]
                  & ((pUf[patchi][facei].x() - U[celli].x())*d);

                gradUy[celli] +=
                    AinvLocal[celli]
                  & ((pUf[patchi][facei].y() - U[celli].y())*d);

                gradUz[celli] +=
                    AinvLocal[celli]
                  & ((pUf[patchi][facei].z() - U[celli].z())*d);
            }

        }
    }

    forAll(gradU, celli)
    {
        gradU[celli] =
            tensor(gradUx[celli], gradUy[celli], gradUz[celli]);
    }

    volTensorField::Boundary& pgradU(gradU.boundaryFieldRef());
    forAll(pgradU, patchi)
    {
        forAll(pgradU[patchi], facei)
        {
            pgradU[patchi][facei] =
                tensor
                (
                    gradUx.boundaryField()[patchi][facei],
                    gradUy.boundaryField()[patchi][facei],
                    gradUz.boundaryField()[patchi][facei]
                );
        }
    }
    // gradU.correctBoundaryConditions();
    return tgradU;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::reconstruct
(
    GeometricField<scalar, fvPatchField, volMesh>& U,
    const volVectorField& gradU,
    GeometricField<scalar, fvsPatchField, surfaceMesh>& UOwn,
    GeometricField<scalar, fvsPatchField, surfaceMesh>& UNei
)
{
    autoPtr<ReconstructionScheme<scalar>> ULimiter
    (
        ReconstructionScheme<scalar>::New(U, U.name(), true)
    );
    UOwn = ULimiter->interpolateOwn();
    UNei = ULimiter->interpolateNei();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::reconstruct
(
    volVectorField& U,
    const volTensorField& gradU,
    surfaceVectorField& UOwn,
    surfaceVectorField& UNei
)
{
    autoPtr<ReconstructionScheme<vector>> ULimiter
    (
        ReconstructionScheme<vector>::New(U, U.name(), true)
    );
    UOwn = ULimiter->interpolateOwn();
    UNei = ULimiter->interpolateNei();
return;

    // const volVectorField& C = mesh_.C();
    // const surfaceVectorField& Cf = mesh_.Cf();
    // const labelList& owner = mesh_.owner();
    // const labelList& neighbour = mesh_.neighbour();
    // volVectorField limiter(calcLimiter(U, gradU));
    // forAll(owner, facei)
    // {
    //     const label& own = owner[facei];
    //     const label& nei = neighbour[facei];
    //
    //     UOwn[facei] =
    //         U[own]
    //       + cmptMultiply(limiter[own], (gradU[own] & (Cf[facei] - C[own])));
    //
    //     UNei[facei] =
    //         U[nei]
    //       + cmptMultiply(limiter[nei], (gradU[nei] & (Cf[facei] - C[nei])));
    // }
    //
    // const volVectorField::Boundary& bU = U.boundaryField();
    // const volVectorField::Boundary& blimiter = limiter.boundaryField();
    // const volTensorField::Boundary& bgradU = gradU.boundaryField();
    // surfaceVectorField::Boundary& bUOwn = UOwn.boundaryFieldRef();
    // surfaceVectorField::Boundary& bUNei = UNei.boundaryFieldRef();
    // forAll(bU, patchi)
    // {
    //     const fvPatch& patch = mesh_.boundary()[patchi];
    //     if (bU[patchi].coupled())
    //     {
    //         const vectorField pdOwn(patch.fvPatch::delta());
    //         const vectorField pdNei(pdOwn - patch.delta());
    //
    //         const vectorField pUOwn(bU[patchi].patchInternalField());
    //         const vectorField pUNei(bU[patchi].patchNeighbourField());
    //
    //         const tensorField pgradUOwn(bgradU[patchi].patchInternalField());
    //         const tensorField pgradUNei(bgradU[patchi].patchNeighbourField());
    //
    //         const vectorField plimOwn(blimiter[patchi].patchInternalField());
    //         const vectorField plimNei(blimiter[patchi].patchNeighbourField());
    //
    //         forAll(pdOwn, facei)
    //         {
    //             bUOwn[patchi][facei] =
    //                 pUOwn[facei]
    //               + cmptMultiply
    //                 (
    //                     plimOwn[facei],
    //                     (pgradUOwn[facei] & pdOwn[facei])
    //                 );
    //             bUNei[patchi][facei] =
    //                 pUNei[facei]
    //               + cmptMultiply
    //                 (
    //                     plimNei[facei],
    //                     (pgradUNei[facei] & pdNei[facei])
    //                 );
    //         }
    //     }
    //     else
    //     {
    //         const vectorField pd(patch.delta());
    //         forAll(pd, facei)
    //         {
    //             const label& celli =
    //                 mesh_.boundaryMesh()[patchi].faceCells()[facei];
    //             bUOwn[patchi][facei] =
    //                 U[celli]
    //               + cmptMultiply
    //                 (
    //                     limiter[celli],
    //                     (gradU[celli] & pd[facei])
    //                 );
    //         }
    //     }
    // }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::reconstruct
(
    volTensorField& U,
    const volTensorField& gradUx,
    const volTensorField& gradUy,
    const volTensorField& gradUz,
    GeometricField<tensor, fvsPatchField, surfaceMesh>& UOwn,
    GeometricField<tensor, fvsPatchField, surfaceMesh>& UNei
)
{
    volVectorField Ux(ops_.decomposeTensorX(U));
    volVectorField Uy(ops_.decomposeTensorY(U));
    volVectorField Uz(ops_.decomposeTensorZ(U));

    surfaceVectorField UxOwn(ops_.decomposeTensorX(UOwn));
    surfaceVectorField UyOwn(ops_.decomposeTensorY(UOwn));
    surfaceVectorField UzOwn(ops_.decomposeTensorZ(UOwn));
    surfaceVectorField UxNei(ops_.decomposeTensorX(UNei));
    surfaceVectorField UyNei(ops_.decomposeTensorY(UNei));
    surfaceVectorField UzNei(ops_.decomposeTensorZ(UNei));

    gradientSchemes::reconstruct(Ux, gradUx, UxOwn, UxNei);
    gradientSchemes::reconstruct(Uy, gradUy, UyOwn, UyNei);
    gradientSchemes::reconstruct(Uz, gradUz, UzOwn, UzNei);

    UOwn.replace(tensor::XX, UxOwn.component(vector::X));
    UOwn.replace(tensor::XY, UxOwn.component(vector::Y));
    UOwn.replace(tensor::XZ, UxOwn.component(vector::Z));

    UOwn.replace(tensor::YX, UyOwn.component(vector::X));
    UOwn.replace(tensor::YY, UyOwn.component(vector::Y));
    UOwn.replace(tensor::YZ, UyOwn.component(vector::Z));

    UOwn.replace(tensor::ZX, UzOwn.component(vector::X));
    UOwn.replace(tensor::ZY, UzOwn.component(vector::Y));
    UOwn.replace(tensor::ZZ, UzOwn.component(vector::Z));

    UNei.replace(tensor::XX, UxNei.component(vector::X));
    UNei.replace(tensor::XY, UxNei.component(vector::Y));
    UNei.replace(tensor::XZ, UxNei.component(vector::Z));

    UNei.replace(tensor::YX, UyNei.component(vector::X));
    UNei.replace(tensor::YY, UyNei.component(vector::Y));
    UNei.replace(tensor::YZ, UyNei.component(vector::Z));

    UNei.replace(tensor::ZX, UzNei.component(vector::X));
    UNei.replace(tensor::ZY, UzNei.component(vector::Y));
    UNei.replace(tensor::ZZ, UzNei.component(vector::Z));

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
