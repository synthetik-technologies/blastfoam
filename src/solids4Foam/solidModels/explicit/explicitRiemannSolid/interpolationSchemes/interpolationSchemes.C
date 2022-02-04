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

#include "interpolationSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(interpolationSchemes, 0);


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

interpolationSchemes::interpolationSchemes(const fvMesh& mesh)
:
    mesh_(mesh),
    own_(mesh_.owner()),
    C_(mesh_.C()),
    Cf_(mesh_.Cf()),
    points_(mesh_.points())
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

interpolationSchemes::~interpolationSchemes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<volVectorField> interpolationSchemes::surfaceToVol
(
    const surfaceVectorField& Uf,
    const pointVectorField& pointU
) const
{
    tmp<volVectorField> tvf_v
    (
        volVectorField::New
        (
            "surfaceToVol(" + Uf.name() + ")",
            mesh_,
            dimensioned<vector>("0", Uf.dimensions(), pTraits<vector>::zero)
        )
    );
    volVectorField& U = tvf_v.ref();

    Field<scalar> w(U.size(), Zero);

    forAll(own_, facei)
    {
        const label own = own_[facei];
        const label nei = mesh_.neighbour()[facei];

        const scalar wOwn(1.0/mag(Cf_[facei] - C_[own]));
        const scalar wNei(1.0/mag(Cf_[facei] - C_[nei]));

        U[own] += Uf[facei]*wOwn;
        U[nei] += Uf[facei]*wNei;

        w[own] += wOwn;
        w[nei] += wNei;
    }

    volVectorField::Boundary& pU = U.boundaryFieldRef();
    forAll(U.boundaryField(), patchi)
    {
        bool fix = pointU.boundaryField()[patchi].fixesValue();
        const fvPatch& patch = mesh_.boundary()[patchi];
        const scalarField pw(1.0/mag(patch.fvPatch::delta()));
        pU[patchi] = Uf.boundaryField()[patchi];
        forAll(U.boundaryField()[patchi], facei)
        {
            const label celli =
                mesh_.boundaryMesh()[patchi].faceCells()[facei];

            U[celli] += Uf.boundaryField()[patchi][facei]*pw[facei];
            w[celli] += pw[facei];

            if (fix)
            {
                const label faceID =
                    mesh_.boundary()[patchi].start() + facei;

                forAll(mesh_.faces()[faceID], nodei)
                {
                    const label nodeID = mesh_.faces()[faceID][nodei];
                    scalar nw(1.0/mag(points_[nodeID] - C_[celli]));
                    U[celli] += pointU[nodeID]*nw;
                    w[celli] += nw;

                    for (label i = 0; i < 7; i++)
                    {
                        scalar si(i);
                        nw =
                            1.0
                           /mag
                            (
                                (
                                    ((si + 1)*points_[nodeID])
                                  + (
                                        (7.0 - si)
                                       *Cf_.boundaryField()[patchi][facei]
                                    )
                                )/8.0 - C_[celli]
                            );

                        U[celli] += pointU[nodeID]*nw;
                        w[celli] += nw;
                    }
                }
            }
        }
    }
    U.primitiveFieldRef() /= w;
    U.correctBoundaryConditions();

    return tvf_v;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void interpolationSchemes::pushUntransformedData
(
    List<Type>& pointData
) const
{
    const globalMeshData& gmd = mesh_.globalData();
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void interpolationSchemes::addSeparated
(
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
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

void interpolationSchemes::volToPoint
(
    const volVectorField& U,
    const volTensorField& gradU,
    pointVectorField& pointU,
    const bool interpBoundaries
) const
{

    tmp<pointVectorField> tpointUOld;
    if (!interpBoundaries)
    {
        tpointUOld = pointVectorField::New("pointUOld", pointU);
    }
    if (Pstream::parRun())
    {
        pointScalarField sum
        (
            IOobject
            (
                "volPointSum",
                U.instance(),
                mesh_
            ),
            pointMesh::New(mesh_),
            dimensionedScalar("zero", dimless, 0.0)
        );
        pointU = Zero;

        forAll (points_, nodei)
        {
            forAll (mesh_.pointCells()[nodei], ci)
            {
                const label celli = mesh_.pointCells()[nodei][ci];
                const vector d = points_[nodei] - C_[celli];
                const vector recons = U[celli] + (gradU[celli] & d);

                pointU[nodei] += recons;
                sum[nodei] += 1.0;
            }
        }

        pointConstraints::syncUntransformedData
        (
            mesh_,
            sum,
            plusEqOp<scalar>()
        );
        addSeparated(sum);
        pushUntransformedData(sum);

        forAll(points_, nodei)
        {
            pointU[nodei] = pointU[nodei] / sum[nodei];
        }

        pointConstraints::syncUntransformedData
        (
            mesh_,
            pointU,
            plusEqOp<vector>()
        );
        addSeparated(pointU);
        pushUntransformedData(pointU);
    }

    else
    {
        forAll (mesh_.pointCells(), nodei)
        {
            vector sum = vector::zero;
            label weights = 0;

            forAll (mesh_.pointCells()[nodei], ci)
            {
                const label celli = mesh_.pointCells()[nodei][ci];
                const vector d = points_[nodei] - C_[celli];
                const vector recons = U[celli] + (gradU[celli] & d);

                sum += recons;
                weights++;
            }

            pointU[nodei] = sum/scalar(weights);
        }
    }
    if (!interpBoundaries)
    {
        const pointVectorField& pointUOld = tpointUOld();
        pointVectorField::Boundary& bpointU = pointU.boundaryFieldRef();
        forAll(mesh_.boundaryMesh(), patchi)
        {
            if (!mesh_.boundaryMesh()[patchi].coupled())
            {
                bpointU[patchi].setInInternalField
                (
                    pointU.primitiveFieldRef(),
                    pointUOld.boundaryField()[patchi].patchInternalField()()
                );
            }
        }
    }
//     else
//     {
//         const faceList& faces = mesh_.faces();
//         forAll(mesh_.boundaryMesh(), patchi)
//         {
//             const polyPatch& patch = mesh_.boundaryMesh()[patchi];
//             vectorField sum(patch.nPoints(), Zero);
//             scalarField weights(patch.nPoints(), 0.0);
//             forAll(patch, fi)
//             {
//                 const label facei = patch.start() + fi;
//                 const labelList& fp = faces[facei];
//                 forAll(fp, ni)
//                 {
//                     const label nodei = fp[ni];
//                     const vector d = points_[nodei] - _[celli];
//                     const vector recons = U[celli] + (gradU[celli] & d);
//
//                     sum += recons;
//                     weights++;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<surfaceVectorField> interpolationSchemes::pointToSurface
(
    const pointVectorField& pointU
) const
{
    vector sum(Zero);
    scalar weights = 0.0;

    tmp<surfaceVectorField > tUf
    (
        surfaceVectorField::New
        (
            "pointToSurface(" + pointU.name() + ")",
            mesh_,
            dimensioned<vector>("0", pointU.dimensions(), Zero)
        )
    );
    surfaceVectorField& Uf = tUf.ref();

    forAll(own_, facei)
    {
        sum = vector::zero;
        weights = 0.0;

        forAll(mesh_.faces()[facei], ni)
        {
            const label pointi = mesh_.faces()[facei][ni];
            scalar w = 1.0;
//             scalar w(1.0/mag(mesh_.points()[pointi] - mesh_.Cf()[facei]));
            sum += pointU[pointi]*w;
            weights += w;
        }

        Uf[facei] = sum/weights;
    }

    surfaceVectorField::Boundary& pUf(Uf.boundaryFieldRef());
    forAll(pUf, patchi)
    {
        forAll(pUf[patchi], fi)
        {
            const label& facei = mesh_.boundary()[patchi].start() + fi;
            sum = vector::zero;
            weights = 0.0;

            forAll(mesh_.faces()[facei], ni)
            {
                const label pointi = mesh_.faces()[facei][ni];
                scalar w = 1.0;
//                 scalar w
//                 (
//                     1.0
//                    /mag
//                     (
//                         mesh_.points()[pointi] - mesh_.faceCentres()[facei]
//                     )
//                 );
                sum += pointU[pointi]*w;
                weights += w;
            }

            pUf[patchi][fi] = sum/weights;
        }
    }

    return tUf;
}


tmp<volVectorField> interpolationSchemes::pointToVol
(
    const pointVectorField& pointU
) const
{
    vector sum(Zero);
    scalar weights = 0.0;

    tmp<volVectorField > tU
    (
        volVectorField::New
        (
            "pointToVol(" + pointU.name() + ")",
            mesh_,
            dimensioned<vector>("0", pointU.dimensions(), Zero)
        )
    );
    volVectorField& U = tU.ref();

    forAll(mesh_.cells(), celli)
    {
        sum = vector::zero;
        weights = 0.0;

        labelList points(mesh_.cells()[celli].labels(mesh_.faces()));

        forAll(points, pI)
        {
            const label pointi = points[pI];
            scalar w = 1.0;
//             scalar w(1.0/mag(mesh_.points()[pointi] - mesh_.C()[celli]));
            sum += pointU[pointi]*w;
            weights += w;
        }

        U[celli] = sum/weights;
    }

    volVectorField::Boundary& pU(U.boundaryFieldRef());
    forAll(pU, patchi)
    {
        forAll(pU[patchi], fi)
        {
            const label& facei = mesh_.boundary()[patchi].start() + fi;
            sum = vector::zero;
            weights = 0.0;

            forAll(mesh_.faces()[facei], ni)
            {
                const label pointi = mesh_.faces()[facei][ni];
                scalar w = 1.0;
//                 scalar w
//                 (
//                     1.0
//                    /mag
//                     (
//                         mesh_.points()[pointi] - mesh_.faceCentres()[facei]
//                     )
//                 );
                sum += pointU[pointi]*w;
                weights += w;
            }

            pU[patchi][fi] = sum/weights;
        }
    }

    return tU;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
