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

interpolationSchemes::interpolationSchemes(const fvMesh& vm)
:
    mesh_(vm),
    own_(mesh_.owner()),
    X_(mesh_.C()),
    XF_(mesh_.Cf()),
    XN_(mesh_.points())
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

interpolationSchemes::~interpolationSchemes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volVectorField> interpolationSchemes::surfaceToVol
(
    const GeometricField<vector, fvsPatchField, surfaceMesh>& sf
) const
{
    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf_v
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "surfaceToVol("+sf.name()+')',
                sf.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensioned<vector>("0", sf.dimensions(), pTraits<vector>::zero)
        )
    );
    GeometricField<vector, fvPatchField, volMesh>& U = tvf_v.ref();

    tmp<GeometricField<scalar, fvPatchField, volMesh> > tvf_s
    (
        new GeometricField<scalar, fvPatchField, volMesh>
        (
            IOobject
            (
                "surfaceToVol("+sf.name()+')',
                sf.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensioned<scalar>("0", dimless, pTraits<scalar>::zero)
        )
    );
    GeometricField<scalar, fvPatchField, volMesh>& w = tvf_s.ref();

    forAll(own_, faceID)
    {
        const label ownID = own_[faceID];
        const label neiID = mesh_.neighbour()[faceID];

        const vector dOwn = XF_[faceID] - X_[ownID];
        const vector dNei = XF_[faceID] - X_[neiID];

        U[ownID] += sf[faceID]*(1.0/mag(dOwn));
        U[neiID] += sf[faceID]*(1.0/mag(dNei));

        w[ownID] += 1.0/mag(dOwn);
        w[neiID] += 1.0/mag(dNei);
    }

    const pointVectorField& pointRhoU =
        mesh_.lookupObject<pointVectorField> ("pointRhoU");
    forAll(mesh_.boundary(), patchID)
    {
        forAll(mesh_.boundary()[patchID], facei)
        {
            const label bCellID =
                mesh_.boundaryMesh()[patchID].faceCells()[facei];

            vector d = XF_.boundaryField()[patchID][facei] - X_[bCellID];

            U[bCellID] += sf.boundaryField()[patchID][facei] * (1.0/mag(d));
            w[bCellID] += 1.0/mag(d);

            if (pointRhoU.boundaryField()[patchID].fixesValue())
            {
                const label faceID =
                    mesh_.boundary()[patchID].start() + facei;

                forAll(mesh_.faces()[faceID], nodei)
                {
                    const label nodeID = mesh_.faces()[faceID][nodei];
                    d = XN_[nodeID] - X_[bCellID];
                    U[bCellID] += pointRhoU[nodeID]*(1.0/mag(d));
                    w[bCellID] += (1.0/mag(d));

                    for (int i=0; i<7; i++)
                    {
                        d =
                            ((((i+1)*XN_[nodeID])
                          + ((7-i)*XF_.boundaryField()[patchID][facei]) )/8.0)
                          - X_[bCellID];

                        U[bCellID] += pointRhoU[nodeID]*(1.0/mag(d));
                        w[bCellID] += 1.0/mag(d);
                    }
                }
            }
        }
    }
    U.primitiveFieldRef() /= w.internalField();

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
    const GeometricField<tensor, fvPatchField, volMesh>& Ugrad,
    GeometricField<vector, pointPatchField, pointMesh>& Un
) const
{

    const fvMesh& mesh = mesh_;
    const volVectorField& C = mesh.C();


    if( Pstream::parRun() )
    {
        pointScalarField sum
        (
            IOobject
            (
                "volPointSum",
                mesh.polyMesh::instance(),
                mesh
            ),
            pointMesh::New(mesh),
            dimensionedScalar("zero", dimless, 0.0)
        );

        forAll (mesh.points(), nodeID)
        {
            Un[nodeID] = vector::zero;
        }

        forAll (mesh.points(), nodeID)
        {
            forAll (mesh.pointCells()[nodeID], cell)
            {
                const label& cellID = mesh.pointCells()[nodeID][cell];
                const vector& d = mesh.points()[nodeID] - C[cellID];
                const vector& recons = U[cellID] + ( Ugrad[cellID] & d );
                const scalar& weight = 1;

                Un[nodeID] += recons;
                sum[nodeID] += weight;
            }
        }

        pointConstraints::syncUntransformedData(mesh, sum, plusEqOp<scalar>());
        addSeparated(sum);
        pushUntransformedData(sum);

        forAll (mesh.points(), nodeID)
        {
            Un[nodeID] = Un[nodeID] / sum[nodeID];
        }

        pointConstraints::syncUntransformedData(mesh, Un, plusEqOp<vector>());
        addSeparated(Un);
        pushUntransformedData(Un);
    }

    else
    {
        forAll (mesh.pointCells(), nodeID)
        {
            vector sum = vector::zero;
            scalar weights = 0.0;

            forAll (mesh.pointCells()[nodeID], cell)
            {
                const label& cellID = mesh.pointCells()[nodeID][cell];
                const vector& d = mesh.points()[nodeID] - C[cellID];
                const vector& recons = U[cellID] + ( Ugrad[cellID] & d );

                //sum += recons * (1.0/mag(d));
                //weights += (1.0/mag(d));

                sum += recons;
                weights += 1.0;
            }

            Un[nodeID] = sum / weights;
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

    forAll(own_, faceID)
    {
        sum = vector::zero;
        weights = 0.0;

        forAll(mesh_.faces()[faceID], node)
        {
            const label& nodeID = mesh_.faces()[faceID][node];
            // d = XN_[nodeID] - XF_[faceID];
            // sum += U[nodeID]*(1.0/mag(d));
            // weights += 1.0/mag(d);
            sum += U[nodeID];
            weights += 1.0;
        }

        Uf[faceID] = sum/weights;
    }

    forAll(mesh_.boundary(), patchID)
    {
        forAll(mesh_.boundary()[patchID], facei)
        {
            const label& faceID = mesh_.boundary()[patchID].start() + facei;
            sum = vector::zero;
            weights = 0.0;

            forAll(mesh_.faces()[faceID], node)
            {
                const label& nodeID = mesh_.faces()[faceID][node];
                // d = XN_[nodeID] - XF_.boundaryField()[patchID][facei];
                // sum += U[nodeID]*(1.0/mag(d));
                // weights += 1.0/mag(d);
                sum += U[nodeID];
                weights += 1.0;
            }

            Uf.boundaryFieldRef()[patchID][facei] = sum/weights;
        }
    }

    return tsf_v;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
