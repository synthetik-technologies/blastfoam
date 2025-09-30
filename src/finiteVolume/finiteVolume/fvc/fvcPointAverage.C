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

#include "fvcPointAverage.H"
#include "valuePointPatchField.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Pure average

// surface to volume averaging see fvcAverage.H

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> average
(
    const GeometricField<Type, pointPatchField, pointMesh>& vpf
)
{
    const fvMesh& mesh = dynamicCast<const fvMesh>(vpf.mesh().mesh());
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tvsf
    (
        GeometricField<Type, fvsPatchField, surfaceMesh>::New
        (
            "average(" + vpf.name() + ")",
            mesh,
            dimensioned<Type>(vpf.dimensions(), Zero)
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& vsf = tvsf.ref();

    const Field<Type>& vpfI = vpf.primitiveField();

    const faceList& faces = mesh.faces();
    forAll(vsf, facei)
    {
        const face& f = faces[facei];
        forAll(f, pi)
        {
            vsf[facei] += vpfI[f[pi]];
        }

        vsf[facei] /= scalar(f.size());
    }

    typename GeometricField<Type, fvsPatchField, surfaceMesh>::Boundary& bvsf =
        vsf.boundaryFieldRef();
    const typename GeometricField<Type, pointPatchField, pointMesh>::Boundary& bvpf =
        vpf.boundaryField();
    forAll(bvsf, patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];
        Field<Type>& pvsf = bvsf[patchi];

        if (bvpf[patchi].fixesValue())
        {
            const Field<Type>& pvpf =
                dynamicCast<const valuePointPatchField<Type>>(bvpf[patchi]);
            forAll(pvsf, fi)
            {
                const face& f = patch.localFaces()[fi];
                forAll(f, pi)
                {
                    pvsf[fi] += pvpf[f[pi]];
                }
                pvsf[fi] /= scalar(f.size());
            }
        }
        else
        {
            forAll(pvsf, fi)
            {
                const face& f = patch[fi];
                forAll(f, pi)
                {
                    pvsf[fi] += vpfI[f[pi]];
                }
                pvsf[fi] /= scalar(f.size());
            }
        }
    }

    return tvsf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> pointVolAverage
(
    const GeometricField<Type, pointPatchField, pointMesh>& vpf
)
{
    const fvMesh& mesh = dynamicCast<const fvMesh>(vpf.mesh().mesh());
    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf
    (
        GeometricField<Type, fvsPatchField, volMesh>::New
        (
            "average(" + vpf.name() + ")",
            mesh,
            dimensioned<Type>(vpf.dimensions(), Zero)
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf.ref();

    const faceList& faces = mesh.faces();
    const labelListList& cellPoints = mesh.cellPoints();
    forAll(cellPoints, celli)
    {
        const labelList& cp = cellPoints[celli];
        forAll(cp, pi)
        {
            vf[celli] += vpf[cp[pi]];
        }

        vf[celli] /= scalar(cp.size());
    }

    typename GeometricField<Type, fvPatchField, surfaceMesh>::Boundary& bvf =
        vf.boundaryFieldRef();
    const typename GeometricField<Type, pointPatchField, pointMesh>::Boundary& bvpf =
        vpf.boundaryField();
    if (Pstream::parRun())
    {
        Field<Type> bValues(mesh.nFaces() - mesh.nInternalFaces());
        Field<scalar> bWeights(mesh.nFaces() - mesh.nInternalFaces(), 1.0);
        forAll(bvf, patchi)
        {
            const polyPatch& patch = mesh.boundaryMesh()[patchi];
            const label start = patch.start() - mesh.nInternalFaces();
            Field<Type>& pvf = bvf[patchi];

            forAll(pvf, fi)
            {
                const face& f = patch[fi];

                forAll(f, pi)
                {
                    pvf[fi] += vpf[f[pi]];
                }

                bValues[start + fi] = pvf[fi];
                bWeights[start + fi] = scalar(f.size());
            }
        }
        syncTools::syncBoundaryFaceList(mesh, bValues, plusEqOp<Type>());
        syncTools::syncBoundaryFaceList(mesh, bWeights, plusEqOp<scalar>());

        forAll(bvf, patchi)
        {
            const polyPatch& patch = mesh.boundaryMesh()[patchi];
            const label start = patch.start() - mesh.nInternalFaces();
            Field<Type>& pvf = bvf[patchi];

            forAll(pvf, fi)
            {
                pvf[fi] = bValues[start + fi]/bWeights[start + fi];
            }
        }
    }
    else
    {
        forAll(bvf, patchi)
        {
            const polyPatch& patch = mesh.boundaryMesh()[patchi];
            Field<Type>& pvf = bvf[patchi];

            forAll(pvf, fi)
            {
                const face& f = patch[fi];
                forAll(f, pi)
                {
                    const label pointi = f[pi];
                    pvf[fi] += vpf[pointi];
                }

                pvf[fi] /= scalar(f.size());
            }
        }
    }
    vf.correctBoundaryConditions();

    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
