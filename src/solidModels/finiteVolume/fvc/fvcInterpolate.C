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

#include "fvcInterpolate.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace pointFieldOps
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Helper functions

template<class Type>
void pushUntransformedData(List<Type>& pointData, const polyMesh& mesh)
{
    const globalMeshData& gmd = mesh.globalData();
    const indirectPrimitivePatch& cpp = gmd.coupledPatch();
    const labelList& meshPoints = cpp.meshPoints();

    const mapDistribute& slavesMap = gmd.globalCoPointSlavesMap();
    const labelListList& slaves = gmd.globalCoPointSlaves();

    List<Type> elems(slavesMap.constructSize());
    forAll(meshPoints, i)
    {
        elems[i] = pointData[meshPoints[i]];
    }

    forAll(slaves, i)
    {
        const labelList& slavePoints = slaves[i];
        forAll(slavePoints, j)
        {
            elems[slavePoints[j]] = elems[i];
        }
    }

    slavesMap.reverseDistribute(elems.size(), elems, false);

    forAll(meshPoints, i)
    {
        pointData[meshPoints[i]] = elems[i];
    }
}


template<class Type>
void addSeparated(GeometricField<Type, pointPatchField, pointMesh>& pf)
{
    typename GeometricField<Type, pointPatchField, pointMesh>::
        Internal& pfi = pf.ref();

    typename GeometricField<Type, pointPatchField, pointMesh>::
        Boundary& pfbf = pf.boundaryFieldRef();

    forAll(pfbf, patchi)
    {
        if (pfbf[patchi].coupled())
        {
            refCast<coupledPointPatchField<Type>>
                (pfbf[patchi]).initSwapAddSeparated
                (
                    Pstream::commsTypes::nonBlocking,
                    pfi
                );
        }
    }

    Pstream::waitRequests();

    forAll(pfbf, patchi)
    {
        if (pfbf[patchi].coupled())
        {
            refCast<coupledPointPatchField<Type>>
                (pfbf[patchi]).swapAddSeparated
                (
                    Pstream::commsTypes::nonBlocking,
                    pfi
                );
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pointFieldOps

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void volPointInterpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    >& vfGrad,
    GeometricField<Type, pointPatchField, pointMesh>& vpf,
    const bool boundary
)
{
    const fvMesh& mesh = vf.mesh();
    const volVectorField& C = mesh.C();
    const pointField& points = mesh.points();
    const labelListList& pointCells = mesh.pointCells();

    Field<Type>& vpfI = vpf.primitiveFieldRef();

    labelHashSet fixedPoints;
    const typename GeometricField<Type, pointPatchField, pointMesh>::Boundary& bvpf =
        vpf.boundaryField();
    forAll(bvpf, patchi)
    {
        if (bvpf[patchi].fixesValue())
        {
            fixedPoints.insert(mesh.boundaryMesh()[patchi].meshPoints());
        }
    }

    if (Pstream::parRun())
    {
        pointScalarField sumWeights
        (
            IOobject
            (
                "volPointSum",
                mesh.polyMesh::instance(),
                mesh
            ),
            pointMesh::New(mesh),
            dimensionedScalar("zero", dimless, Zero)
        );

        forAll(points, pointi)
        {
            vpfI[pointi] = Zero;
            const labelList& pc = pointCells[pointi];
            forAll(pc, ci)
            {
                const label celli = pc[ci];
                const vector d(points[pointi] - C[celli]);
                const scalar w = 1.0/mag(d);

                vpfI[pointi] += (vf[celli] + (vfGrad[celli] & d))*w;
                sumWeights[pointi] += w;
            }
        }

        if (boundary)
        {
            forAll(mesh.boundary(), patchi)
            {
                const polyPatch& patch = mesh.boundaryMesh()[patchi];
                const vectorField& pC = C.boundaryField()[patchi];
                const fvPatchField<Type>& pvf = vf.boundaryField()[patchi];
                const fvPatchField<typename outerProduct<Type, vector>::type>&
                    pvfGrad = vfGrad.boundaryField()[patchi];
                forAll(pvf, fi)
                {
                    const face& f = patch[fi];
                    forAll(f, pi)
                    {
                        const label pointi = f[pi];
                        if (!fixedPoints.found(pointi))
                        {
                            const vector d(points[pointi] - pC[fi]);
                            const scalar w = 1.0/mag(d);

                            vpfI[pointi] += (pvf[fi] + (pvfGrad[fi] & d))*w;
                            sumWeights[pointi] += w;
                        }
                    }
                }
            }
        }

        pointConstraints::syncUntransformedData(mesh, sumWeights, plusEqOp<scalar>());
        pointFieldOps::addSeparated(sumWeights);
        pointFieldOps::pushUntransformedData(sumWeights, mesh);

        vpfI /= sumWeights;

        pointConstraints::syncUntransformedData(mesh, vpf, plusEqOp<vector>());
        pointFieldOps::addSeparated(vpf);
        pointFieldOps::pushUntransformedData(vpf, mesh);
    }

    else
    {
        scalarField sumWeights(vpf.size(), Zero);
        forAll(points, pointi)
        {
            vpfI[pointi] = Zero;
            const labelList& pc = pointCells[pointi];
            forAll(pc, ci)
            {
                const label celli = pc[ci];
                const vector d(points[pointi] - C[celli]);
                const scalar w = 1.0/mag(d);

                vpfI[pointi] += (vf[celli] + (vfGrad[celli] & d))*w;
                sumWeights[pointi] += w;
            }
        }
        if (boundary)
        {
            forAll(mesh.boundary(), patchi)
            {
                if (!vpf.boundaryField()[patchi].fixesValue())
                {
                    const polyPatch& patch = mesh.boundaryMesh()[patchi];
                    const vectorField& pC = C.boundaryField()[patchi];
                    const fvPatchField<Type>& pvf = vf.boundaryField()[patchi];
                    const fvPatchField<typename outerProduct<Type, vector>::type>&
                        pvfGrad = vfGrad.boundaryField()[patchi];
                    forAll(patch, fi)
                    {
                        const face& f = patch[fi];
                        forAll(f, pi)
                        {
                            const label pointi = f[pi];
                            const vector d(points[pointi] - pC[fi]);
                            const scalar w = 1.0/mag(d);

                            vpfI[pointi] += (pvf[fi] + (pvfGrad[fi] & d))*w;
                            sumWeights[pointi] += w;
                        }
                    }
                }
            }
        }
        vpfI /= sumWeights;
    }
    vpf.correctBoundaryConditions();
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> surfVolInterpolate
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& vsf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
    (
        GeometricField<Type, fvPatchField, volMesh>::New
        (
            "surfaceToVol(" + vsf.name() + ")",
            vsf.mesh(),
            dimensioned<Type>(vsf.dimensions(), Zero)
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf.ref();

    Field<scalar> weights(vf.size(), Zero);

    const fvMesh& mesh = vsf.mesh();
    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();
    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();
    forAll(vsf, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        const scalar wOwn(1.0/mag(Cf[facei] - C[own]));
        const scalar wNei(1.0/mag(Cf[facei] - C[nei]));

        vf[own] += vsf[facei]*wOwn;
        vf[nei] += vsf[facei]*wNei;

        weights[own] += wOwn;
        weights[nei] += wNei;
    }

    typename GeometricField<Type, fvPatchField, volMesh>::Boundary& bvf =
        vf.boundaryFieldRef();
    const surfaceVectorField::Boundary& bvsf = vsf.boundaryField();
    forAll(bvf, patchi)
    {
        const fvPatch& patch = vsf.mesh().boundary()[patchi];
        const scalarField pw(patch.fvPatch::deltaCoeffs());

        forAll(bvf[patchi], facei)
        {
            const label celli = patch.faceCells()[facei];

            vf[celli] += bvsf[patchi][facei]*pw[facei];
            weights[celli] += pw[facei];
        }

        // Set boundary field
        bvf[patchi] = bvsf[patchi];
    }

    vf.primitiveFieldRef() /= weights;
    vf.correctBoundaryConditions();

    return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> pointSurfInterpolate
(
    const GeometricField<Type, pointPatchField, pointMesh>& vpf
)
{
    const fvMesh& mesh = dynamicCast<const fvMesh>(vpf.mesh().mesh());
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tvsf
    (
        GeometricField<Type, fvsPatchField, surfaceMesh>::New
        (
            "pointToSurface(" + vpf.name() + ")",
            mesh,
            dimensioned<Type>(vpf.dimensions(), Zero)
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& vsf = tvsf.ref();

    const surfaceVectorField& Cf = mesh.Cf();
    const pointField& points = mesh.points();
    const faceList& faces = mesh.faces();
    forAll(vsf, facei)
    {
        scalar sumW = 0.0;
        const face& f = faces[facei];
        const vector& cf = Cf[facei];

        forAll(f, pi)
        {
            const label pointi = f[pi];
            scalar w(1.0/mag(points[pointi] - cf));
            vsf[facei] += vpf[pointi]*w;
            sumW += w;
        }

        vsf[facei] /= sumW;
    }

    typename GeometricField<Type, fvsPatchField, surfaceMesh>::Boundary& bvsf =
        vsf.boundaryFieldRef();
    if (Pstream::parRun())
    {
        Field<Type> bValues(mesh.nFaces() - mesh.nInternalFaces());
        Field<scalar> bWeights(mesh.nFaces() - mesh.nInternalFaces(), 1.0);
        forAll(bvsf, patchi)
        {
            const polyPatch& patch = mesh.boundaryMesh()[patchi];
            const label start = patch.start() - mesh.nInternalFaces();
            fvsPatchField<Type>& pvsf = bvsf[patchi];
            const vectorField& pCf = Cf.boundaryField()[patchi];

            forAll(pCf, fi)
            {
                const face& f = patch[fi];
                scalar sumW = 0.0;

                forAll(f, pi)
                {
                    const label pointi = f[pi];
                    scalar w(1.0/mag(points[pointi] - pCf[fi]));
                    pvsf[fi] += vpf[pointi]*w;
                    sumW += w;
                }

                bValues[start + fi] = pvsf[fi];
                bWeights[start + fi] = sumW;
            }
        }
        syncTools::syncBoundaryFaceList(mesh, bValues, plusEqOp<Type>());
        syncTools::syncBoundaryFaceList(mesh, bWeights, plusEqOp<scalar>());

        forAll(bvsf, patchi)
        {
            const polyPatch& patch = mesh.boundaryMesh()[patchi];
            const label start = patch.start() - mesh.nInternalFaces();
            fvsPatchField<Type>& pvsf = bvsf[patchi];

            forAll(pvsf, fi)
            {
                pvsf[fi] = bValues[start + fi]/bWeights[start + fi];
            }
        }
    }
    else
    {
        forAll(bvsf, patchi)
        {
            const polyPatch& patch = mesh.boundaryMesh()[patchi];
            fvsPatchField<Type>& pvsf = bvsf[patchi];
            const vectorField& pCf = Cf.boundaryField()[patchi];

            forAll(pCf, fi)
            {
                const face& f = patch[fi];
                scalar sumW = 0.0;

                forAll(f, pi)
                {
                    const label pointi = f[pi];
                    scalar w(1.0/mag(points[pointi] - pCf[fi]));
                    pvsf[fi] += vpf[pointi]*w;
                    sumW += w;
                }

                pvsf[fi] /= sumW;
            }
        }
    }

    return tvsf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> pointVolInterpolate
(
    const GeometricField<Type, pointPatchField, pointMesh>& vpf
)
{
    const fvMesh& mesh = dynamicCast<const fvMesh>(vpf.mesh().mesh());
    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf
    (
        GeometricField<Type, fvsPatchField, volMesh>::New
        (
            "pointToVol(" + vpf.name() + ")",
            mesh,
            dimensioned<Type>(vpf.dimensions(), Zero)
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf.ref();

    const volVectorField& C = mesh.C();
    const pointField& points = mesh.points();
    const labelListList& cellPoints = mesh.cellPoints();
    forAll(vf, celli)
    {
        const labelList& cp = cellPoints[celli];
        scalar sumW = 0.0;
        forAll(cp, pi)
        {
            const label pointi = cp[pi];
            scalar w(1.0/mag(points[pointi] - C[celli]));
            vf[celli] += vpf[pointi]*w;
            sumW += w;
        }

        vf[celli] /= sumW;
    }

    typename GeometricField<Type, fvPatchField, volMesh>::Boundary& bvf =
        vf.boundaryFieldRef();
    if (Pstream::parRun())
    {
        Field<Type> bValues(mesh.nFaces() - mesh.nInternalFaces());
        Field<scalar> bWeights(mesh.nFaces() - mesh.nInternalFaces(), 1.0);
        forAll(bvf, patchi)
        {
            const polyPatch& patch = mesh.boundaryMesh()[patchi];
            const label start = patch.start() - mesh.nInternalFaces();
            Field<Type>& pvf = bvf[patchi];
            const vectorField& pC = C.boundaryField()[patchi];

            forAll(pvf, fi)
            {
                const face& f = patch[fi];
                scalar sumW = 0.0;

                forAll(f, pi)
                {
                    const label pointi = f[pi];
                    scalar w(1.0/mag(points[pointi] - pC[fi]));
                    pvf[fi] += vpf[pointi]*w;
                    sumW += w;
                }

                bValues[start + fi] = pvf[fi];
                bWeights[start + fi] = sumW;
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
            const vectorField& pC = C.boundaryField()[patchi];

            forAll(pvf, fi)
            {
                const face& f = patch[fi];
                scalar sumW = 0.0;

                forAll(f, pi)
                {
                    const label pointi = f[pi];
                    scalar w(1.0/mag(points[pointi] - pC[fi]));
                    pvf[fi] += vpf[pointi]*w;
                    sumW += w;
                }

                pvf[fi] /= sumW;
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
