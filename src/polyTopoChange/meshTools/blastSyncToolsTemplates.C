/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "blastSyncTools.H"
#include "fvPatchField.H"
#include "valuePointPatchField.H"
#include "surfaceFields.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class GeoField>
void Foam::blastSyncTools::correctProcessorBoundaries(const fvMesh& mesh)
{
    HashTable<GeoField*> flds
    (
        const_cast<fvMesh&>(mesh).lookupClass<GeoField>()
    );
    forAllIter(typename HashTable<GeoField*>, flds, iter)
    {
        GeoField& fld = *iter();

        //mimic "evaluate" but only for coupled patches (processor or cyclic)
        // and only for blocking or nonBlocking comms (no scheduled comms)
        if
        (
            Pstream::defaultCommsType == Pstream::commsTypes::blocking
         || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
        )
        {
            label nReq = Pstream::nRequests();

            forAll(fld.boundaryField(), patchi)
            {
                if (isA<processorPolyPatch>(mesh.boundaryMesh()[patchi]))
                {
                    fld.boundaryFieldRef()[patchi].initEvaluate
                    (
                        Pstream::defaultCommsType
                    );
                }
            }

            // Block for any outstanding requests
            if
            (
                Pstream::parRun()
             && Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
            )
            {
                Pstream::waitRequests(nReq);
            }

            forAll(fld.boundaryField(), patchi)
            {
                if (isA<processorPolyPatch>(mesh.boundaryMesh()[patchi]))
                {
                    fld.boundaryFieldRef()[patchi].evaluate
                    (
                        Pstream::defaultCommsType
                    );
                }
            }
        }
        else
        {
            //Scheduled patch updates not supported
            FatalErrorInFunction
                << "Unsuported communications type "
                << Pstream::commsTypeNames[Pstream::defaultCommsType]
                << exit(FatalError);
        }
    }
}


template<class Type>
void Foam::blastSyncTools::setInPointBoundaries(const fvMesh& mesh)
{
    typedef GeometricField<Type, pointPatchField, pointMesh> GeoField;
    HashTable<GeoField*> flds
    (
        const_cast<fvMesh&>(mesh).lookupClass<GeoField>()
    );
    forAllIter(typename HashTable<GeoField*>, flds, iter)
    {
        GeoField& fld = *iter();
        typename GeoField::Boundary& bfld = fld.boundaryFieldRef();
        forAll(bfld, patchi)
        {
            if (isA<valuePointPatchField<Type>>(bfld[patchi]))
            {
                bfld[patchi].setInternalField
                (
                    fld.primitiveFieldRef(),
                    dynamicCast<const Field<Type>>(bfld[patchi])
                );
            }
        }
    }
}


template<class Type>
void Foam::blastSyncTools::correctPointBoundaries(const fvMesh& mesh)
{
    typedef GeometricField<Type, pointPatchField, pointMesh> GeoField;
    HashTable<GeoField*> flds
    (
        const_cast<fvMesh&>(mesh).lookupClass<GeoField>()
    );
    forAllIter(typename HashTable<GeoField*>, flds, iter)
    {
        GeoField& fld = *iter();
        typename GeoField::Boundary& bfld = fld.boundaryFieldRef();
        forAll(bfld, patchi)
        {
            if (isA<valuePointPatchField<Type>>(bfld[patchi]))
            {
                bfld[patchi] == bfld[patchi].patchInternalField();
            }
        }
    }
}


template<class Type>
void Foam::blastSyncTools::pushUntransformedData
(
    const polyMesh& mesh,
    Field<Type>& pointData
)
{
    const globalMeshData& gmd = mesh.globalData();
//     Field<scalar> nSharedPoints(pointData.size(), 1);
//     gmd.syncPointData
//     (
//         pointData,
//         plusEqOp<Type>(),
//         distributionMap::transform()
//     );
//     gmd.syncPointData
//     (
//         nSharedPoints,
//         plusEqOp<scalar>(),
//         distributionMap::transform()
//     );
//     pointData /= nSharedPoints;

    // Transfer onto coupled patch

    const indirectPrimitivePatch& cpp = gmd.coupledPatch();
    const labelList& meshPoints = cpp.meshPoints();

    const distributionMap& slavesMap = gmd.globalCoPointSlavesMap();
    const labelListList& slaves = gmd.globalCoPointSlaves();

    List<Type> elems(slavesMap.constructSize());
    forAll(meshPoints, i)
    {
        elems[i] = pointData[meshPoints[i]];
    }

    // Combine master data with slave data
    forAll(slaves, i)
    {
        const labelList& slavePoints = slaves[i];

        // Copy master data to slave slots
        forAll(slavePoints, j)
        {
            elems[slavePoints[j]] = elems[i];
        }
    }

    // Push slave-slot data back to slaves
    slavesMap.reverseDistribute(elems.size(), elems, false);

    // Extract back onto mesh
    forAll(meshPoints, i)
    {
        pointData[meshPoints[i]] = elems[i];
    }
}

// ************************************************************************* //
