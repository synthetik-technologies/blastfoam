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
    own_(mesh_.owner()),
    nei_(mesh_.neighbour()),
    C_(mesh_.C()),
    Cf_(mesh_.Cf()),
    points_(mesh_.points()),

    Ainv_
    (
        IOobject("Ainv", mesh_),
        mesh_,
        dimensionedTensor("Ainv", dimensionSet(0,2,0,0,0,0,0), tensor::zero)
    ),

    AinvLocal_
    (
        IOobject("AinvLocal", mesh_),
        mesh_,
        dimensionedTensor
        (
            "AinvLocal",
            dimensionSet(0,2,0,0,0,0,0),
            tensor::zero
        )
    )
{
    Ainv_ = distanceMatrix();
    AinvLocal_ = distanceMatrixLocal();
}


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

tmp<volTensorField> gradientSchemes::distanceMatrix() const
{
    tmp<GeometricField<tensor, fvPatchField, volMesh> > tAinv
    (
        GeometricField<tensor, fvPatchField, volMesh>::New
        (
            "Ainv",
            mesh_,
            dimensioned<tensor>("0", sqr(dimLength), pTraits<tensor>::zero)
        )
    );
    volTensorField& Ainv = tAinv.ref();

    forAll(own_, facei)
    {
        const label own = own_[facei];
        const label nei = nei_[facei];
        const tensor dSqr = (C_[nei] - C_[own])*(C_[nei] - C_[own]);

        Ainv[own] += dSqr;
        Ainv[nei] += dSqr;
    }

    forAll(mesh_.boundary(), patchi)
    {
        if (mesh_.boundary()[patchi].coupled())
        {
            vectorField pCNei
            (
                C_.boundaryField()[patchi].patchNeighbourField()
            );

            forAll(mesh_.boundary()[patchi], facei)
            {
                const label celli =
                    mesh_.boundaryMesh()[patchi].faceCells()[facei];

                const vector d = pCNei[facei] - C_[celli];
                Ainv[celli] += d*d;
            }
        }
    }

    if (Pstream::parRun())
    {
        Ainv.correctBoundaryConditions();
    }

    Ainv.primitiveFieldRef() = inv(Ainv.primitiveField());

    return tAinv;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volTensorField> gradientSchemes::distanceMatrixLocal() const
{
    const volVectorField& rhoU = mesh_.lookupObject<volVectorField> ("rhoU");

    tmp<GeometricField<tensor, fvPatchField, volMesh> > tAinv
    (
        GeometricField<tensor, fvPatchField, volMesh>::New
        (
            "Ainv",
            mesh_,
            dimensioned<tensor>("0", sqr(dimLength), pTraits<tensor>::zero)
        )
    );
    GeometricField<tensor, fvPatchField, volMesh>& Ainv = tAinv.ref();

    tmp<GeometricField<tensor, fvPatchField, volMesh> > tdCd
    (
        GeometricField<tensor, fvPatchField, volMesh>::New
        (
            "distanceMatrixLocal",
            mesh_,
            dimensioned<tensor>("0", sqr(dimLength), pTraits<tensor>::zero)
        )
    );
    GeometricField<tensor, fvPatchField, volMesh>& dCd = tdCd.ref();

    forAll(own_, facei)
    {
        const label own = own_[facei];
        const label nei = nei_[facei];
        const vector dOwn = Cf_[facei] - C_[own];
        const vector dNei = Cf_[facei] - C_[nei];

        dCd[own] += dOwn*dOwn;
        dCd[nei] += dNei*dNei;
    }

    forAll(mesh_.boundary(), patchi)
    {
        if (rhoU.boundaryField()[patchi].fixesValue())
        {
            forAll(mesh_.boundary()[patchi], facei)
            {
                const label celli =
                    mesh_.boundaryMesh()[patchi].faceCells()[facei];
                const label mfacei =
                    mesh_.boundary()[patchi].start() + facei;

                vector d = Cf_.boundaryField()[patchi][facei] - C_[celli];
                dCd[celli] += d*d;

                forAll(mesh_.faces()[mfacei], nodei)
                {
                    const label pti = mesh_.faces()[mfacei][nodei];

                    d = points_[pti] - C_[celli];
                    dCd[celli] += d*d;

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
                        dCd[celli] += d*d;
                    }
                }
            }
        }
    }

    Ainv.primitiveFieldRef() = inv(dCd.primitiveField());
    return tAinv;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volVectorField> gradientSchemes::gradient
(
    const GeometricField<scalar, fvPatchField, volMesh>& U
)   const
{
    tmp<GeometricField<vector, fvPatchField, volMesh> > tgradU
    (
        GeometricField<vector, fvPatchField, volMesh>::New
        (
            "gradient("+U.name()+')',
            mesh_,
            dimensioned<vector>
            (
                "0",
                U.dimensions()/dimLength,
                pTraits<vector>::zero
            )
        )
    );
    GeometricField<vector, fvPatchField, volMesh>& gradU = tgradU.ref();

    forAll(own_, faceID)
    {
        const label own = own_[faceID];
        const label nei = nei_[faceID];
        const vector d = C_[nei] - C_[own];

        gradU[own] += Ainv_[own] & (U[nei] - U[own])*d;
        gradU[nei] -= Ainv_[nei] & (U[own] - U[nei])*d;
    }


    forAll(mesh_.boundary(), patchi)
    {
        if (mesh_.boundary()[patchi].coupled())
        {
            vectorField CNei
            (
                C_.boundaryField()[patchi].patchNeighbourField()
            );

            scalarField UNei
            (
                U.boundaryField()[patchi].patchNeighbourField()
            );

            forAll(mesh_.boundary()[patchi], facei)
            {
                const label& celli =
                    mesh_.boundaryMesh()[patchi].faceCells()[facei];

                const vector d = CNei[facei] - C_[celli];

                gradU[celli] += Ainv_[celli] & (UNei[facei] - U[celli])*d;
            }
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
    const GeometricField<vector, fvPatchField, volMesh>& U
)   const
{
    GeometricField<vector, fvPatchField, volMesh> gradUx
    (
        gradientSchemes::gradient(U.component(0))
    );
    GeometricField<vector, fvPatchField, volMesh> gradUy
    (
        gradientSchemes::gradient(U.component(1))
    );
    GeometricField<vector, fvPatchField, volMesh> gradUz
    (
        gradientSchemes::gradient(U.component(2))
    );

    tmp<GeometricField<tensor, fvPatchField, volMesh> > tgradU
    (
        GeometricField<tensor, fvPatchField, volMesh>::New
        (
            "gradient(" + U.name() + ')',
            mesh_,
            dimensioned<tensor>
            (
                "0",
                U.dimensions()/dimLength,
                pTraits<tensor>::zero
            )
        )
    );
    GeometricField<tensor, fvPatchField, volMesh>& gradU = tgradU.ref();

    forAll(gradU, celli)
    {
        gradU[celli] = tensor(gradUx[celli], gradUy[celli], gradUz[celli]);
    }

    if( Pstream::parRun() )
    {
        gradU.correctBoundaryConditions();
    }

    return tgradU;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::gradient
(
    const GeometricField<tensor, fvPatchField, volMesh>& U,
    GeometricField<tensor, fvPatchField, volMesh>& UgradX,
    GeometricField<tensor, fvPatchField, volMesh>& UgradY,
    GeometricField<tensor, fvPatchField, volMesh>& UgradZ
)   const
{
    GeometricField<vector, fvPatchField, volMesh> Ux
    (
        GeometricField<vector, fvPatchField, volMesh>::New
        (
            U.name() + 'x',
            mesh_,
            dimensioned<vector>("0", U.dimensions(), pTraits<vector>::zero)
        )
    );
    GeometricField<vector, fvPatchField, volMesh> Uy
    (
        GeometricField<vector, fvPatchField, volMesh>::New
        (
            U.name() + 'y',
            mesh_,
            dimensioned<vector>("0", U.dimensions(), pTraits<vector>::zero)
        )
    );
    GeometricField<vector, fvPatchField, volMesh> Uz
    (
        GeometricField<vector, fvPatchField, volMesh>::New
        (
            U.name() + 'z',
            mesh_,
            dimensioned<vector>("0", U.dimensions(), pTraits<vector>::zero)
        )
    );

    operations(mesh_).decomposeTensor(U, Ux, Uy, Uz);

    if (Pstream::parRun())
    {
        Ux.correctBoundaryConditions();
        Uy.correctBoundaryConditions();
        Uz.correctBoundaryConditions();
    }

    UgradX = gradientSchemes::gradient(Ux);
    UgradY = gradientSchemes::gradient(Uy);
    UgradZ = gradientSchemes::gradient(Uz);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volTensorField> gradientSchemes::localGradient
(
    const GeometricField<vector, fvPatchField, volMesh>& U,
    const GeometricField<vector, fvsPatchField, surfaceMesh>& Uf
) const
{
    const objectRegistry& db = mesh_.thisDb();
    const pointVectorField& pointRhoU = db.lookupObject<pointVectorField> ("pointRhoU");

    tmp<GeometricField<tensor, fvPatchField, volMesh> > tgradU
    (
        GeometricField<tensor, fvPatchField, volMesh>::New
        (
            "gradient("+U.name()+')',
            mesh_,
            dimensioned<tensor>
            (
                "0",
                U.dimensions()/dimLength,
                Zero
            )
        )
    );
    GeometricField<tensor, fvPatchField, volMesh>& gradU = tgradU.ref();

    vectorField gradUx(gradU.size(), Zero);
    vectorField gradUy(gradU.size(), Zero);
    vectorField gradUz(gradU.size(), Zero);

    forAll(own_, facei)
    {
        const label own = own_[facei];
        const label nei = nei_[facei];
        const vector dOwn = Cf_[facei] - C_[own];
        const vector dNei = Cf_[facei] - C_[nei];

        gradUx[own] += AinvLocal_[own] & ((Uf[facei].x() - U[own].x())*dOwn);
        gradUy[own] += AinvLocal_[own] & ((Uf[facei].y() - U[own].y())*dOwn);
        gradUz[own] += AinvLocal_[own] & ((Uf[facei].z() - U[own].z())*dOwn);

        gradUx[nei] += AinvLocal_[nei] & ((Uf[facei].x() - U[nei].x())*dNei);
        gradUy[nei] += AinvLocal_[nei] & ((Uf[facei].y() - U[nei].y())*dNei);
        gradUz[nei] += AinvLocal_[nei] & ((Uf[facei].z() - U[nei].z())*dNei);
    }

    forAll(mesh_.boundary(), patchi)
    {
        forAll(mesh_.boundary()[patchi], facei)
        {
            const label celli =
                mesh_.boundaryMesh()[patchi].faceCells()[facei];

            vector d = Cf_.boundaryField()[patchi][facei] - C_[celli];

            gradUx[celli] +=
                AinvLocal_[celli]
              & ((Uf.boundaryField()[patchi][facei].x() - U[celli].x())*d);

            gradUy[celli] +=
                AinvLocal_[celli]
              & ((Uf.boundaryField()[patchi][facei].y() - U[celli].y())*d);

            gradUz[celli] +=
                AinvLocal_[celli]
              & ((Uf.boundaryField()[patchi][facei].z() - U[celli].z())*d);

            if (pointRhoU.boundaryField()[patchi].fixesValue())
            {
                const label mfacei =
                    mesh_.boundary()[patchi].start() + facei;

                forAll(mesh_.faces()[mfacei], pti)
                {
                    const label mpti = mesh_.faces()[mfacei][pti];
                    vector d = points_[mpti] - C_[celli];

                    gradUx[celli] +=
                        AinvLocal_[celli]
                      & ((pointRhoU[mpti].x() - U[celli].x())*d);

                    gradUy[celli] +=
                        AinvLocal_[celli]
                      & ((pointRhoU[mpti].y() - U[celli].y())*d);

                    gradUz[celli] +=
                        AinvLocal_[celli]
                      & ((pointRhoU[mpti].z() - U[celli].z())*d);

                    for (label i = 0; i < 7; i++)
                    {
                        scalar si(i);
                        d =
                            (
                                (
                                    (si + 1.0)*points_[mpti]
                                  + (
                                        (7.0 - si)
                                       *Cf_.boundaryField()[patchi][facei]
                                    )
                                )/8.0
                            ) - C_[celli];

                        gradUx[celli] +=
                            AinvLocal_[celli]
                          & ((pointRhoU[mpti].x() - U[celli].x())*d);

                        gradUy[celli] +=
                            AinvLocal_[celli]
                          & ((pointRhoU[mpti].y() - U[celli].y())*d);

                        gradUz[celli] +=
                            AinvLocal_[celli]
                          & ((pointRhoU[mpti].z() - U[celli].z())*d);
                    }
                }
            }
        }
    }

    forAll(gradU, celli)
    {
        gradU[celli] = tensor(gradUx[celli], gradUy[celli], gradUz[celli]);
    }

    return tgradU;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::reconstruct
(
    GeometricField<scalar, fvPatchField, volMesh>& U,
    const GeometricField<vector, fvPatchField, volMesh>& gradU,
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
        forAll(mesh_.boundaryMesh()[patchi], facei)
        {
            const label& celli =
                mesh_.boundaryMesh()[patchi].faceCells()[facei];

            U.boundaryFieldRef()[patchi][facei] =
                U[celli]
              + (
                    gradU[celli]
                  & (Cf_.boundaryField()[patchi][facei] - C_[celli])
                );

            UOwn.boundaryFieldRef()[patchi][facei] =
                U[celli]
              + (
                    gradU[celli]
                  & (Cf_.boundaryField()[patchi][facei] - C_[celli])
                );
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::reconstruct
(
    GeometricField<vector, fvPatchField, volMesh>& U,
    const GeometricField<tensor, fvPatchField, volMesh>& gradU,
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
        forAll(mesh_.boundaryMesh()[patchi], facei)
        {
            const label& celli =
                mesh_.boundaryMesh()[patchi].faceCells()[facei];

            U.boundaryFieldRef()[patchi][facei] =
                U[celli]
              + (
                    gradU[celli]
                  & (Cf_.boundaryField()[patchi][facei] - C_[celli])
                );
            UOwn.boundaryFieldRef()[patchi][facei] =
                U[celli]
              + (
                    gradU[celli]
                  & (Cf_.boundaryField()[patchi][facei] - C_[celli])
                );
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::reconstruct
(
    GeometricField<tensor, fvPatchField, volMesh>& U,
    const GeometricField<tensor, fvPatchField, volMesh>& gradUx,
    const GeometricField<tensor, fvPatchField, volMesh>& gradUy,
    const GeometricField<tensor, fvPatchField, volMesh>& gradUz,
    GeometricField<tensor, fvsPatchField, surfaceMesh>& UOwn,
    GeometricField<tensor, fvsPatchField, surfaceMesh>& UNei
)
{
    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf
    (
        GeometricField<vector, fvPatchField, volMesh>::New
        (
            "reconstruct("+U.name()+')',
            mesh_,
            dimensionedVector(U.dimensions(), Zero)
        )
    );
    GeometricField<vector, fvPatchField, volMesh> Ux(tvf());
    GeometricField<vector, fvPatchField, volMesh> Uy(tvf());
    GeometricField<vector, fvPatchField, volMesh>& Uz(tvf.ref());

    operations op(mesh_);
    op.decomposeTensor(U, Ux, Uy, Uz);

    tmp<GeometricField<vector, fvsPatchField, surfaceMesh> > tUxOwn
    (
        GeometricField<vector, fvsPatchField, surfaceMesh>::New
        (
            "reconstruct(" + U.name() + ')',
            mesh_,
            dimensionedVector(U.dimensions(), Zero)
        )
    );
    GeometricField<vector, fvsPatchField, surfaceMesh>& UxOwn(tUxOwn.ref());
    GeometricField<vector, fvsPatchField, surfaceMesh> UyOwn(UxOwn);
    GeometricField<vector, fvsPatchField, surfaceMesh> UzOwn(UxOwn);
    GeometricField<vector, fvsPatchField, surfaceMesh> UxNei(UxOwn);
    GeometricField<vector, fvsPatchField, surfaceMesh> UyNei(UxOwn);
    GeometricField<vector, fvsPatchField, surfaceMesh> UzNei(UxOwn);

//     op.decomposeTensor(UOwn, UxOwn, UyOwn, UzOwn);
//     op.decomposeTensor(UNei, UxNei, UyNei, UzNei);

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
        forAll(mesh_.boundaryMesh()[patchi], facei)
        {
            const label& celli =
                mesh_.boundaryMesh()[patchi].faceCells()[facei];

            UOwn.boundaryFieldRef()[patchi][facei] =
                tensor
                (
                    Ux[celli]
                  + (
                        gradUx[celli]
                      & (Cf_.boundaryField()[patchi][facei] - C_[celli])
                    ),
                    Uy[celli]
                  + (
                        gradUy[celli]
                      & (Cf_.boundaryField()[patchi][facei] - C_[celli])
                    ),
                    Uz[celli]
                  + (
                        gradUz[celli]
                      & (Cf_.boundaryField()[patchi][facei] - C_[celli])
                    )
                );
        }
        U.boundaryFieldRef()[patchi] == UOwn.boundaryField()[patchi];
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
