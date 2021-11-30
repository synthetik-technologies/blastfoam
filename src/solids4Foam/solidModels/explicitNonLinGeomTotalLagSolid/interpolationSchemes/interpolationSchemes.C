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
    const GeometricField<vector, fvsPatchField, surfaceMesh>& Uf,
    const GeometricField<vector, pointPatchField, pointMesh>& pU
) const
{
    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf_v
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "surfaceToVol("+Uf.name()+')',
                Uf.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensioned<vector>("0", Uf.dimensions(), pTraits<vector>::zero)
        )
    );
    GeometricField<vector, fvPatchField, volMesh>& U = tvf_v.ref();

    Field<scalar> w(U.size(), Zero);

    forAll(own_, facei)
    {
        const label own = own_[facei];
        const label nei = mesh_.neighbour()[facei];

        const vector dOwn = Cf_[facei] - C_[own];
        const vector dNei = Cf_[facei] - C_[nei];

        U[own] += Uf[facei]*(1.0/mag(dOwn));
        U[nei] += Uf[facei]*(1.0/mag(dNei));

        w[own] += 1.0/mag(dOwn);
        w[nei] += 1.0/mag(dNei);
    }

    forAll(U.boundaryField(), patchi)
    {
        forAll(U.boundaryField()[patchi], facei)
        {
            const label celli =
                mesh_.boundaryMesh()[patchi].faceCells()[facei];

            vector d = Cf_.boundaryField()[patchi][facei] - C_[celli];

            U[celli] += Uf.boundaryField()[patchi][facei]/mag(d);
            w[celli] += 1.0/mag(d);

            if (pU.boundaryField()[patchi].fixesValue())
            {
                const label faceID =
                    mesh_.boundary()[patchi].start() + facei;

                forAll(mesh_.faces()[faceID], nodei)
                {
                    const label nodeID = mesh_.faces()[faceID][nodei];
                    d = points_[nodeID] - C_[celli];
                    U[celli] += pU[nodeID]/mag(d);
                    w[celli] += (1.0/mag(d));

                    for (label i = 0; i < 7; i++)
                    {
                        scalar si(i);
                        d =
                            (
                                (
                                    ((si + 1)*points_[nodeID])
                                  + (
                                        (7.0 - si)
                                       *Cf_.boundaryField()[patchi][facei]
                                    )
                                )/8.0
                            ) - C_[celli];

                        U[celli] += pU[nodeID]/mag(d);
                        w[celli] += 1.0/mag(d);
                    }
                }
            }
        }
    }
    U.primitiveFieldRef() /= w;

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
    const GeometricField<vector, fvPatchField, volMesh>& U,
    const GeometricField<tensor, fvPatchField, volMesh>& gradU,
    GeometricField<vector, pointPatchField, pointMesh>& pU
) const
{

    if( Pstream::parRun() )
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

        pU.ref() = Zero;
        forAll (points_, nodeID)
        {
            forAll (mesh_.pointCells()[nodeID], cell)
            {
                const label cellID = mesh_.pointCells()[nodeID][cell];
                const vector d = points_[nodeID] - C_[cellID];
                const vector recons = U[cellID] + (gradU[cellID] & d);

                pU[nodeID] += recons;
                sum[nodeID] += 1.0;
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

        forAll(points_, nodeID)
        {
            pU[nodeID] = pU[nodeID] / sum[nodeID];
        }

        pointConstraints::syncUntransformedData
        (
            mesh_,
            pU,
            plusEqOp<vector>()
        );
        addSeparated(pU);
        pushUntransformedData(pU);
    }

    else
    {
        forAll (mesh_.pointCells(), nodeID)
        {
            vector sum = vector::zero;
            scalar weights = 0.0;

            forAll (mesh_.pointCells()[nodeID], cell)
            {
                const label cellID = mesh_.pointCells()[nodeID][cell];
                const vector d = points_[nodeID] - C_[cellID];
                const vector recons = U[cellID] + (gradU[cellID] & d);

                //sum += recons * (1.0/mag(d));
                //weights += (1.0/mag(d));

                sum += recons;
                weights += 1.0;
            }

            pU[nodeID] = sum/weights;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<surfaceVectorField> interpolationSchemes::pointToSurface
(
    const GeometricField<vector, pointPatchField, pointMesh>& U
) const
{
    // vector d = vector::zero;
    vector sum = vector::zero;
    scalar weights = 0.0;

    tmp<GeometricField<vector, fvsPatchField, surfaceMesh> > tsf_v
    (
        new GeometricField<vector, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "pointToSurface("+U.name()+')',
                U.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensioned<vector>("0", U.dimensions(), pTraits<vector>::zero)
        )
    );
    GeometricField<vector, fvsPatchField, surfaceMesh>& Uf = tsf_v.ref();

    forAll(own_, facei)
    {
        sum = vector::zero;
        weights = 0.0;

        forAll(mesh_.faces()[facei], node)
        {
            const label& nodeID = mesh_.faces()[facei][node];
            // d = XN_[nodeID] - XF_[facei];
            // sum += U[nodeID]*(1.0/mag(d));
            // weights += 1.0/mag(d);
            sum += U[nodeID];
            weights += 1.0;
        }

        Uf[facei] = sum/weights;
    }

    forAll(mesh_.boundary(), patchi)
    {
        forAll(mesh_.boundary()[patchi], facei)
        {
            const label& faceID = mesh_.boundary()[patchi].start() + facei;
            sum = vector::zero;
            weights = 0.0;

            forAll(mesh_.faces()[faceID], node)
            {
                const label& nodeID = mesh_.faces()[faceID][node];
                // d = XN_[nodeID] - XF_.boundaryField()[patchi][facei];
                // sum += U[nodeID]*(1.0/mag(d));
                // weights += 1.0/mag(d);
                sum += U[nodeID];
                weights += 1.0;
            }

            Uf.boundaryFieldRef()[patchi][facei] = sum/weights;
        }
    }

    return tsf_v;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
