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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(gradientSchemes, 0);


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

gradientSchemes::gradientSchemes
(
    const fvMesh& vm
)
:
    MeshObject<fvMesh, MoveableMeshObject, gradientSchemes>(vm),
    mesh_(vm),
    ops_(mesh_),
    own_(mesh_.owner()),
    nei_(mesh_.neighbour()),
    C_(mesh_.C()),
    Cf_(mesh_.Cf()),
    points_(mesh_.points()),

    Ainv_(distanceMatrix()),
    AinvLocal_(distanceMatrixLocal())
{}


bool gradientSchemes::movePoints()
{
    Ainv_ = distanceMatrix();
    AinvLocal_ = distanceMatrixLocal();
    return true;
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

gradientSchemes::~gradientSchemes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tensor gradientSchemes::stabInv(const tensor& t)
{
    if (magSqr(t) < small)
    {
        return tensor::zero;
    }

    scalar scale = magSqr(t);
    Vector<bool> removeCmpts
    (
        magSqr(t.xx())/scale < small,
        magSqr(t.yy())/scale < small,
        magSqr(t.zz())/scale < small
    );
    if (removeCmpts.x() || removeCmpts.y() || removeCmpts.z())
    {
        tensor tPlus(t);

        if (removeCmpts.x())
        {
            tPlus += tensor(1,0,0,0,0,0,0,0,0);
        }

        if (removeCmpts.y())
        {
            tPlus += tensor(0,0,0,0,1,0,0,0,0);
        }

        if (removeCmpts.z())
        {
            tPlus += tensor(0,0,0,0,0,0,0,0,1);
        }

        tensor tInv = inv(tPlus);

        if (removeCmpts.x())
        {
            tInv -= tensor(1,0,0,0,0,0,0,0,0);
        }

        if (removeCmpts.y())
        {
            tInv -= tensor(0,0,0,0,1,0,0,0,0);
        }

        if (removeCmpts.z())
        {
            tInv -= tensor(0,0,0,0,0,0,0,0,1);
        }
        return tInv;
    }
    return inv(t);
}


tmp<tensorField> gradientSchemes::distanceMatrix() const
{
    tmp<tensorField> tAinv(new tensorField(mesh_.nCells(), Zero));
    tensorField& Ainv = tAinv.ref();

    forAll(own_, facei)
    {
        const label own = own_[facei];
        const label nei = nei_[facei];
        const vector dOwn = C_[nei] - C_[own];
        const vector dNei  = C_[own] - C_[nei];

        Ainv[own] += dOwn*dOwn;
        Ainv[nei] += dNei*dNei;
    }

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];
        const vectorField pd(patch.delta());

        if (patch.coupled())
        {
            forAll(mesh_.boundary()[patchi], facei)
            {
                const label celli =
                    mesh_.boundaryMesh()[patchi].faceCells()[facei];
                Ainv[celli] += pd[facei]*pd[facei];
            }
        }
    }

    forAll(Ainv, celli)
    {
        Ainv[celli] = stabInv(Ainv[celli]);
    }

    return tAinv;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<tensorField> gradientSchemes::distanceMatrixLocal() const
{
    tmp<tensorField > tAinv(new tensorField(mesh_.nCells(), Zero));
    tensorField& Ainv = tAinv.ref();

    forAll(own_, facei)
    {
        const label own = own_[facei];
        const label nei = nei_[facei];
        const vector dOwn = Cf_[facei] - C_[own];
        const vector dNei = Cf_[facei] - C_[nei];

        Ainv[own] += dOwn*dOwn;
        Ainv[nei] += dNei*dNei;
    }

    const volVectorField& pD = mesh_.lookupObject<volVectorField> ("D");
    forAll(mesh_.boundary(), patchi)
    {
        bool fix = pD.boundaryField()[patchi].fixesValue();
        const fvPatch& patch = mesh_.boundary()[patchi];
        const vectorField pd(patch.fvPatch::delta());

        forAll(mesh_.boundary()[patchi], facei)
        {
            const label celli =
                mesh_.boundaryMesh()[patchi].faceCells()[facei];

            Ainv[celli] += pd[facei]*pd[facei];

            if (fix)
            {
                const label mfacei =
                    mesh_.boundary()[patchi].start() + facei;
                forAll(mesh_.faces()[mfacei], nodei)
                {
                    const label pti = mesh_.faces()[mfacei][nodei];

                    vector d = points_[pti] - C_[celli];
                    Ainv[celli] += d*d;

                    for (label i = 0; i < 7; i++)
                    {
                        scalar si(i);
                        d =
                            (
                                (
                                    (si + 1.0)*points_[pti]
                                  + (
                                        (7.0 - si)
                                       *Cf_.boundaryField()[patchi][facei]
                                    )
                                )
                            )/8.0 - C_[celli];
                        Ainv[celli] += d*d;
                    }
                }
            }
        }
    }

    return tAinv;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volVectorField> gradientSchemes::gradient
(
    const GeometricField<scalar, fvPatchField, volMesh>& U
)   const
{
    tmp<volVectorField> tgradU
    (
        volVectorField::New
        (
            "grad(" + U.name() + ")",
            mesh_,
            dimensioned<vector>
            (
                "0",
                U.dimensions()/dimLength,
                pTraits<vector>::zero
            )
        )
    );
    volVectorField& gradU = tgradU.ref();

    forAll(own_, faceID)
    {
        const label own = own_[faceID];
        const label nei = nei_[faceID];
        const vector dOwn = C_[nei] - C_[own];
        const vector dNei = C_[own] - C_[nei];

        gradU[own] += Ainv_[own] & ((U[nei] - U[own])*dOwn);
        gradU[nei] += Ainv_[nei] & ((U[own] - U[nei])*dNei);
    }

    const volScalarField::Boundary& pU(U.boundaryField());
    volVectorField::Boundary& pgradU = gradU.boundaryFieldRef();
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];
        if (mesh_.boundary()[patchi].coupled())
        {
            const vectorField pd(mesh_.boundary()[patchi].delta());
            const scalarField UNei
            (
                pU[patchi].patchNeighbourField()
            );
            pgradU[patchi] =
                patch.deltaCoeffs()
               *patch.nf()
               *(UNei - pU[patchi].patchInternalField());


            forAll(mesh_.boundary()[patchi], facei)
            {
                const label& celli =
                    mesh_.boundaryMesh()[patchi].faceCells()[facei];

                gradU[celli] +=
                    Ainv_[celli] & ((UNei[facei] - U[celli])*pd[facei]);
            }
        }
        else
        {
            pgradU[patchi] =
                patch.deltaCoeffs()
               *patch.nf()
               *(pU[patchi] - pU[patchi].patchInternalField());
        }
    }

    if (Pstream::parRun())
    {
        gradU.correctBoundaryConditions();
    }

    return tgradU;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volTensorField> gradientSchemes::gradient
(
    const volVectorField& U
)   const
{
    volVectorField gradUx(gradientSchemes::gradient(U.component(0)));
    volVectorField gradUy(gradientSchemes::gradient(U.component(1)));
    volVectorField gradUz(gradientSchemes::gradient(U.component(2)));

    tmp<volTensorField > tgradU
    (
        volTensorField::New
        (
            "grad(" + U.name() + ")",
            mesh_,
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
    const GeometricField<vector, pointPatchField, pointMesh>& pointU
) const
{
    tmp<volTensorField > tgradU
    (
        volTensorField::New
        (
            "grad(" + U.name() + ")",
            mesh_,
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
            "gradUx",
            mesh_,
            dimensionedVector(gradU.dimensions(), Zero)
        )
    );
    volVectorField gradUy
    (
        volVectorField::New
        (
            "gradUy",
            mesh_,
            dimensionedVector(gradU.dimensions(), Zero)
        )
    );
    volVectorField gradUz
    (
        volVectorField::New
        (
            "gradUz",
            mesh_,
            dimensionedVector(gradU.dimensions(), Zero)
        )
    );

    forAll(own_, facei)
    {
        const label own = own_[facei];
        const label nei = nei_[facei];
        const vector dOwn = Cf_[facei] - C_[own];
        const vector dNei = Cf_[facei] - C_[nei];

        gradUx[own] +=
            AinvLocal_[own] & ((Uf[facei].x() - U[own].x())*dOwn);
        gradUy[own] +=
            AinvLocal_[own] & ((Uf[facei].y() - U[own].y())*dOwn);
        gradUz[own] +=
            AinvLocal_[own] & ((Uf[facei].z() - U[own].z())*dOwn);

        gradUx[nei] +=
            AinvLocal_[nei] & ((Uf[facei].x() - U[nei].x())*dNei);
        gradUy[nei] +=
            AinvLocal_[nei] & ((Uf[facei].y() - U[nei].y())*dNei);
        gradUz[nei] +=
            AinvLocal_[nei] & ((Uf[facei].z() - U[nei].z())*dNei);
    }

    const volVectorField::Boundary& pU(U.boundaryField());
    const surfaceVectorField::Boundary& pUf(Uf.boundaryField());
    volVectorField::Boundary& pgradUx(gradUx.boundaryFieldRef());
    volVectorField::Boundary& pgradUy(gradUy.boundaryFieldRef());
    volVectorField::Boundary& pgradUz(gradUz.boundaryFieldRef());

    forAll(U.boundaryField(), patchi)
    {
        bool fix = pointU.boundaryField()[patchi].fixesValue();
        const fvPatch& patch = mesh_.boundary()[patchi];
        tensorField pgradU
        (
            patch.deltaCoeffs()
           *patch.nf()
           *(pUf[patchi] - pU[patchi].patchInternalField())
        );
        ops_.decomposeTensor
        (
            pgradU,
            pgradUx[patchi],
            pgradUy[patchi],
            pgradUz[patchi]
        );

        forAll(C_.boundaryField()[patchi], facei)
        {
            const label celli =
                mesh_.boundaryMesh()[patchi].faceCells()[facei];
            vector d = Cf_.boundaryField()[patchi][facei] - C_[celli];

            gradUx[celli] +=
                AinvLocal_[celli]
              & ((pUf[patchi][facei].x() - U[celli].x())*d);

            gradUy[celli] +=
                AinvLocal_[celli]
              & ((pUf[patchi][facei].y() - U[celli].y())*d);

            gradUz[celli] +=
                AinvLocal_[celli]
              & ((pUf[patchi][facei].z() - U[celli].z())*d);


            if (fix)
            {
                const label fi =
                    mesh_.boundary()[patchi].start() + facei;

                forAll(mesh_.faces()[fi], pti)
                {
                    const label pi = mesh_.faces()[fi][pti];
                    vector d = points_[pi] - C_[celli];

                    gradUx[celli] +=
                        AinvLocal_[celli]
                      & ((pointU[pi].x() - U[celli].x())*d);

                    gradUy[celli] +=
                        AinvLocal_[celli]
                      & ((pointU[pi].y() - U[celli].y())*d);

                    gradUz[celli] +=
                        AinvLocal_[celli]
                      & ((pointU[pi].z() - U[celli].z())*d);

                    for (label i = 0; i < 7; i++)
                    {
                        scalar si(i);
                        d =
                            (
                                (
                                    (si + 1.0)*points_[pi]
                                  + (
                                        (7.0 - si)
                                       *Cf_.boundaryField()[patchi][facei]
                                    )
                                )/8.0
                            ) - C_[celli];

                        gradUx[celli] +=
                            AinvLocal_[celli]
                          & ((pointU[pi].x() - U[celli].x())*d);

                        gradUy[celli] +=
                            AinvLocal_[celli]
                          & ((pointU[pi].y() - U[celli].y())*d);

                        gradUz[celli] +=
                            AinvLocal_[celli]
                          & ((pointU[pi].z() - U[celli].z())*d);
                    }
                }
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
    if (Pstream::parRun())
    {
        gradU.correctBoundaryConditions();
    }
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
    forAll(own_, facei)
    {
        const label& own = own_[facei];
        const label& nei = nei_[facei];

        UOwn[facei] = U[own] + (gradU[own] & (Cf_[facei] - C_[own]));
        UNei[facei] = U[nei] + (gradU[nei] & (Cf_[facei] - C_[nei]));
    }

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];
        if (patch.coupled())
        {
            const vectorField pdOwn(patch.fvPatch::delta());
            const vectorField pdNei(pdOwn - patch.delta());

            const scalarField pUOwn
            (
                U.boundaryField()[patchi].patchInternalField()
            );
            const scalarField pUNei
            (
                U.boundaryField()[patchi].patchNeighbourField()
            );

            const vectorField pgradUOwn
            (
                gradU.boundaryField()[patchi].patchInternalField()
            );
            const vectorField pgradUNei
            (
                gradU.boundaryField()[patchi].patchNeighbourField()
            );

            forAll(pdOwn, facei)
            {
                UOwn.boundaryFieldRef()[patchi][facei] =
                    pUOwn[facei] + (pgradUOwn[facei] & pdOwn[facei]);
                UNei.boundaryFieldRef()[patchi][facei] =
                    pUNei[facei] + (pgradUNei[facei] & pdNei[facei]);
            }
        }
        else
        {
            const vectorField pd(patch.delta());
            forAll(pd, facei)
            {
                const label& celli =
                    mesh_.boundaryMesh()[patchi].faceCells()[facei];

                UOwn.boundaryFieldRef()[patchi][facei] =
                    U[celli] + (gradU[celli] & pd[facei]);
            }
            UNei.boundaryFieldRef()[patchi] = U.boundaryField()[patchi];
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::reconstruct
(
    volVectorField& U,
    const volTensorField& gradU,
    GeometricField<vector, fvsPatchField, surfaceMesh>& UOwn,
    GeometricField<vector, fvsPatchField, surfaceMesh>& UNei
)
{
    forAll(own_, facei)
    {
        const label& own = own_[facei];
        const label& nei = nei_[facei];

        UOwn[facei] = U[own] + (gradU[own] & (Cf_[facei] - C_[own]));
        UNei[facei] = U[nei] + (gradU[nei] & (Cf_[facei] - C_[nei]));
    }

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];
        if (patch.coupled())
        {
            const vectorField pdOwn(patch.fvPatch::delta());
            const vectorField pdNei(pdOwn - patch.delta());

            const vectorField pUOwn
            (
                U.boundaryField()[patchi].patchInternalField()
            );
            const vectorField pUNei
            (
                U.boundaryField()[patchi].patchNeighbourField()
            );

            const tensorField pgradUOwn
            (
                gradU.boundaryField()[patchi].patchInternalField()
            );
            const tensorField pgradUNei
            (
                gradU.boundaryField()[patchi].patchNeighbourField()
            );

            forAll(pdOwn, facei)
            {
                UOwn.boundaryFieldRef()[patchi][facei] =
                    pUOwn[facei] + (pgradUOwn[facei] & pdOwn[facei]);
                UNei.boundaryFieldRef()[patchi][facei] =
                    pUNei[facei] + (pgradUNei[facei] & pdNei[facei]);
            }
        }
        else
        {
            const vectorField pd(patch.delta());
            forAll(pd, facei)
            {
                const label& celli =
                    mesh_.boundaryMesh()[patchi].faceCells()[facei];

                UOwn.boundaryFieldRef()[patchi][facei] =
                    U[celli] + (gradU[celli] & pd[facei]);
            }
            UNei.boundaryFieldRef()[patchi] = U.boundaryField()[patchi];
        }
    }
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

    forAll(own_, facei)
    {
        UOwn[facei] = tensor(UxOwn[facei], UyOwn[facei], UzOwn[facei]);
        UNei[facei] = tensor(UxNei[facei], UyNei[facei], UzNei[facei]);
    }

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];
        if (patch.coupled())
        {
            const vectorField pdOwn(patch.fvPatch::delta());
            const vectorField pdNei(pdOwn - patch.delta());
            const vectorField pUxOwn
            (
                Ux.boundaryField()[patchi].patchInternalField()
            );
            const vectorField pUxNei
            (
                Ux.boundaryField()[patchi].patchNeighbourField()
            );
            const tensorField pgradUxOwn
            (
                gradUx.boundaryField()[patchi].patchInternalField()
            );
            const tensorField pgradUxNei
            (
                gradUx.boundaryField()[patchi].patchNeighbourField()
            );
            const vectorField pUyOwn
            (
                Uy.boundaryField()[patchi].patchInternalField()
            );
            const vectorField pUyNei
            (
                Uy.boundaryField()[patchi].patchNeighbourField()
            );
            const tensorField pgradUyOwn
            (
                gradUy.boundaryField()[patchi].patchInternalField()
            );
            const tensorField pgradUyNei
            (
                gradUy.boundaryField()[patchi].patchNeighbourField()
            );
            const vectorField pUzOwn
            (
                Uz.boundaryField()[patchi].patchInternalField()
            );
            const vectorField pUzNei
            (
                Uz.boundaryField()[patchi].patchNeighbourField()
            );
            const tensorField pgradUzOwn
            (
                gradUz.boundaryField()[patchi].patchInternalField()
            );
            const tensorField pgradUzNei
            (
                gradUz.boundaryField()[patchi].patchNeighbourField()
            );

            forAll(pdOwn, facei)
            {
                UOwn.boundaryFieldRef()[patchi][facei] =
                    tensor
                    (
                        pUxOwn[facei]
                      + (pgradUxOwn[facei] & pdOwn[facei]),
                        pUyOwn[facei]
                      + (pgradUyOwn[facei] & pdOwn[facei]),
                        pUzOwn[facei]
                      + (pgradUzOwn[facei] & pdOwn[facei])
                    );
                UNei.boundaryFieldRef()[patchi][facei] =
                    tensor
                    (
                        pUxNei[facei]
                      + (pgradUxNei[facei] & pdNei[facei]),
                        pUyNei[facei]
                      + (pgradUyNei[facei] & pdNei[facei]),
                        pUzNei[facei]
                      + (pgradUzNei[facei] & pdNei[facei])
                    );
            }
        }
        else
        {
            const vectorField pd(patch.delta());
            forAll(pd, facei)
            {
                const label& celli =
                    mesh_.boundaryMesh()[patchi].faceCells()[facei];

                UOwn.boundaryFieldRef()[patchi][facei] =
                    tensor
                    (
                        Ux[celli] + (gradUx[celli] & pd[facei]),
                        Uy[celli] + (gradUy[celli] & pd[facei]),
                        Uz[celli] + (gradUz[celli] & pd[facei])
                    );
            }
            UNei.boundaryFieldRef()[patchi] = U.boundaryField()[patchi];
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
