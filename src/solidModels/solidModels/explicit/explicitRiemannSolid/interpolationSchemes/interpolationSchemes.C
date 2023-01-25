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
    mesh_(mesh)
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
    tmp<volVectorField> tU
    (
        volVectorField::New
        (
            "surfaceToVol(" + Uf.name() + ")",
            mesh_,
            dimensioned<vector>(Uf.dimensions(), Zero)
        )
    );
    volVectorField& U = tU.ref();

    Field<scalar> w(U.size(), Zero);

    const labelList& owner = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();
    const volVectorField& C = mesh_.C();
    const surfaceVectorField& Cf = mesh_.Cf();
    const pointField& points = mesh_.points();
    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        const scalar wOwn(1.0/mag(Cf[facei] - C[own]));
        const scalar wNei(1.0/mag(Cf[facei] - C[nei]));

        U[own] += Uf[facei]*wOwn;
        U[nei] += Uf[facei]*wNei;

        w[own] += wOwn;
        w[nei] += wNei;
    }

    volVectorField::Boundary& bU = U.boundaryFieldRef();
    forAll(bU, patchi)
    {
        const pointPatchVectorField& ppointU =
            pointU.boundaryField()[patchi];

        const fvPatch& patch = mesh_.boundary()[patchi];
        const scalarField& pw(patch.deltaCoeffs());
        // bU[patchi] == Uf.boundaryField()[patchi];
        forAll(bU[patchi], facei)
        {
            const label celli =
                mesh_.boundaryMesh()[patchi].faceCells()[facei];

            U[celli] += Uf.boundaryField()[patchi][facei]*pw[facei];
            w[celli] += pw[facei];

            if (ppointU.fixesValue())
            {
                const label faceID =
                    mesh_.boundary()[patchi].start() + facei;

                forAll(mesh_.faces()[faceID], nodei)
                {
                    const label nodeID = mesh_.faces()[faceID][nodei];
                    scalar nw(1.0/mag(points[nodeID] - C[celli]));
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
                                    ((si + 1)*points[nodeID])
                                  + (
                                        (7.0 - si)
                                       *Cf.boundaryField()[patchi][facei]
                                    )
                                )/8.0 - C[celli]
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

    return tU;
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

    const volVectorField& C = mesh_.C();
    const pointField& points = mesh_.points();
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

        forAll(points, nodei)
        {
            forAll(mesh_.pointCells()[nodei], ci)
            {
                const label celli = mesh_.pointCells()[nodei][ci];
                const vector d = points[nodei] - C[celli];
                const vector recons = U[celli] + (gradU[celli] & d);

                // scalar w = 1.0/mag(d);
                scalar w = 1.0;

                pointU[nodei] += recons*w;
                sum[nodei] += w;
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

        forAll(points, nodei)
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
        forAll(mesh_.pointCells(), nodei)
        {
            vector sum = vector::zero;
            label weights = 0;

            forAll(mesh_.pointCells()[nodei], ci)
            {
                const label celli = mesh_.pointCells()[nodei][ci];
                const vector d = points[nodei] - C[celli];
                const vector recons = U[celli] + (gradU[celli] & d);

                // scalar w = 1.0/mag(d);
                scalar w = 1.0;

                sum += recons*w;
                weights += w;
            }

            pointU[nodei] = sum/scalar(weights);
        }
    }
    // if (!interpBoundaries)
    // {
    //     const pointVectorField& pointUOld = tpointUOld();
    //     pointVectorField::Boundary& bpointU = pointU.boundaryFieldRef();
    //     forAll(mesh_.boundaryMesh(), patchi)
    //     {
    //         if (!mesh_.boundaryMesh()[patchi].coupled())
    //         {
    //             bpointU[patchi].setInInternalField
    //             (
    //                 pointU.primitiveFieldRef(),
    //                 pointUOld.boundaryField()[patchi].patchInternalField()()
    //             );
    //         }
    //     }
    // }
    // else
    // {
    //     pointConstraints::setPatchFields(pointU);
    // }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<surfaceVectorField> interpolationSchemes::pointToSurface
(
    const pointVectorField& pointU
) const
{
    tmp<surfaceVectorField > tUf
    (
        surfaceVectorField::New
        (
            "pointToSurface(" + pointU.name() + ")",
            mesh_,
            dimensioned<vector>(pointU.dimensions(), Zero)
        )
    );
    surfaceVectorField& Uf = tUf.ref();

    const surfaceVectorField& Cf = mesh_.Cf();
    forAll(Cf, facei)
    {
        vector sum(Zero);
        scalar weights = 0.0;
        const face& f = mesh_.faces()[facei];
        // const vector& cf = Cf[facei];

        forAll(f, ni)
        {
            const label pointi = f[ni];
            // scalar w(1.0/mag(mesh_.points()[pointi] - cf));
            scalar w = 1.0;
            sum += pointU[pointi]*w;
            weights += w;
        }

        Uf[facei] = sum/weights;
    }

    surfaceVectorField::Boundary& bUf(Uf.boundaryFieldRef());
    forAll(bUf, patchi)
    {
        forAll(bUf[patchi], fi)
        {
            const label& facei = mesh_.boundary()[patchi].start() + fi;
            const face& f = mesh_.faces()[facei];
            // const vector& cf = Cf[facei];
            vector sum(Zero);
            scalar weights = 0.0;

            forAll(f, ni)
            {
                const label pointi = f[ni];
                // scalar w(1.0/mag(mesh_.points()[pointi] - cf));
                scalar w = 1.0;
                sum += pointU[pointi]*w;
                weights += w;
            }

            bUf[patchi][fi] = sum/weights;
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

        const labelList& points(mesh_.cellPoints()[celli]);

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
