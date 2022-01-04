/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

#include "solidSubMeshes.H"
#include "processorFvsPatchField.H"
#include "processorPointPatchFields.H"
// #include "componentMixedPointPatchFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::solidSubMeshes::mapSubMeshVolFields
(
    const PtrList<GeometricField<Type, fvPatchField, volMesh> >& subMeshFields,
    GeometricField<Type, fvPatchField, volMesh>& baseMeshField
) const
{
    // Map internal field from sub-meshes to base mesh

    forAll(subMeshes(), meshI)
    {
        const labelList& cellMap = subMeshes()[meshI].cellMap();

        const GeometricField<Type, fvPatchField, volMesh>& subMeshField =
            subMeshFields[meshI];

        forAll(subMeshField, cellI)
        {
            baseMeshField[cellMap[cellI]] = subMeshField[cellI];
        }

        // Map boundary field

        typename GeometricField<Type, fvPatchField, volMesh>::Boundary&
            pbaseMeshField = baseMeshField.boundaryFieldRef();
        const labelList& patchMap = subMeshes()[meshI].patchMap();
        const labelList& faceMap = subMeshes()[meshI].faceMap();

        forAll(subMeshField.boundaryField(), patchI)
        {
            const fvPatchField<Type>& subMeshFieldP =
                subMeshField.boundaryField()[patchI];

            const label start = subMeshFieldP.patch().patch().start();

            if (patchMap[patchI] != -1)
            {
                fvPatchField<Type>& baseMeshFieldP =
                    pbaseMeshField[patchMap[patchI]];

                if (!baseMeshFieldP.coupled())
                {
                    forAll(subMeshFieldP, faceI)
                    {
                        const label globalGlobalMeshFace =
                            faceMap[start + faceI];

                        const label curGlobalMeshPatchFace =
                            globalGlobalMeshFace
                          - baseMesh().boundaryMesh()[patchMap[patchI]].start();

                        baseMeshFieldP[curGlobalMeshPatchFace] =
                            subMeshFieldP[faceI];
                    }
                }
            }
        }
    }

    baseMeshField.correctBoundaryConditions();
}


template<class Type>
void Foam::solidSubMeshes::mapSubMeshSurfaceFields
(
    const PtrList<GeometricField<Type, fvsPatchField, surfaceMesh> >&
        subMeshFields,
    GeometricField<Type, fvsPatchField, surfaceMesh>& baseMeshField
) const
{
    // Map internal fields from sub-mesh to base mesh

    // Reset field to zero as we will average later
    baseMeshField =
        dimensioned<Type>
        (
            "zero", baseMeshField.dimensions(), pTraits<Type>::zero
        );

    typename GeometricField<Type, fvsPatchField, surfaceMesh>::Boundary&
            pbaseMeshField = baseMeshField.boundaryFieldRef();

    forAll(subMeshes(), meshI)
    {
        const labelList& faceMap = subMeshes()[meshI].faceMap();

        const GeometricField<Type, fvsPatchField, surfaceMesh>& subMeshField =
            subMeshFields[meshI];

        forAll(subMeshField, faceI)
        {
            baseMeshField[faceMap[faceI]] = subMeshField[faceI];
        }

        // Map boundary field

        const labelList& patchMap = subMeshes()[meshI].patchMap();

        forAll(subMeshField.boundaryField(), patchI)
        {
            const fvsPatchField<Type>& subMeshFieldP =
                subMeshField.boundaryField()[patchI];

            const label start = subMeshFieldP.patch().patch().start();

            if (patchMap[patchI] != -1)
            {
                fvsPatchField<Type>& baseMeshFieldP =
                    pbaseMeshField[patchMap[patchI]];

                // Note: unlike volFields, we do map on the coupled patches for
                // surface fields
                forAll(subMeshFieldP, faceI)
                {
                    const label globalGlobalMeshFace = faceMap[start + faceI];

                    const label curGlobalMeshPatchFace =
                        globalGlobalMeshFace
                        - baseMesh().boundaryMesh()[patchMap[patchI]].start();

                    baseMeshFieldP[curGlobalMeshPatchFace] =
                        subMeshFieldP[faceI];
                }
            }
            else // interface faces shared by two materials
            {
                forAll(subMeshFieldP, faceI)
                {
                    const label globalGlobalMeshFace = faceMap[start + faceI];

                    if (globalGlobalMeshFace < baseMesh().nInternalFaces())
                    {
                        // Face value will be the average from both sides
                        baseMeshField[globalGlobalMeshFace] +=
                            0.5*subMeshFieldP[faceI];
                    }
                    else
                    {
                        const label curPatch =
                            baseMesh().boundaryMesh().whichPatch
                            (
                                globalGlobalMeshFace
                            );

                        const label curPatchFace =
                            globalGlobalMeshFace
                            - baseMesh().boundaryMesh()[curPatch].start();

                        pbaseMeshField[curPatch][curPatchFace] =
                            subMeshFieldP[faceI];
                    }
                }
            }
        }
    }

//     baseMeshField.correctBoundaryConditions();

    // Make sure the field is consistent across processor patches
    forAll(baseMeshField.boundaryField(), patchI)
    {
        const fvsPatchField<Type>& baseMeshFieldP =
            baseMeshField.boundaryField()[patchI];

        if (baseMeshFieldP.type() == processorFvsPatchField<Type>::typeName)
        {
            const Field<Type>& patchField = baseMeshFieldP;

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>
                (
                    baseMesh().boundaryMesh()[patchI]
                );

            OPstream::write
            (
                Pstream::commsTypes::blocking,
                procPatch.neighbProcNo(),
                reinterpret_cast<const char*>(patchField.begin()),
                patchField.byteSize()
            );
        }
    }

    forAll(baseMeshField.boundaryField(), patchI)
    {
        fvsPatchField<Type>& baseMeshFieldP =
            pbaseMeshField[patchI];

        if (baseMeshFieldP.type() == processorFvsPatchField<Type>::typeName)
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>
                (
                    baseMesh().boundaryMesh()[patchI]
                );

            Field<Type> ngbPatchField(procPatch.size(), pTraits<Type>::zero);

            IPstream::read
            (
                Pstream::commsTypes::blocking,
                procPatch.neighbProcNo(),
                reinterpret_cast<char*>(ngbPatchField.begin()),
                ngbPatchField.byteSize()
            );

            Field<Type>& patchField = baseMeshFieldP;

            patchField = 0.5*(patchField + ngbPatchField);
        }
    }
}


template<class Type>
void Foam::solidSubMeshes::mapSubMeshPointFields
(
    const PtrList<GeometricField<Type, pointPatchField, pointMesh> >&
        subMeshFields,
    GeometricField<Type, pointPatchField, pointMesh>& baseMeshField
) const
{
    // Map internal field from sub-mesh to base mesh

    vectorField& baseMeshFieldI = baseMeshField.primitiveFieldRef();

    // Reset the baseMeshField to zero as we take global averages later
    baseMeshFieldI = vector::zero;

    // Number of material adjacent to each point
    const labelList& noMat = pointNumOfMaterials();

    // Global point addressing
    const labelList& spLabels =
        baseMesh().globalData().sharedPointLabels();
    const labelList& spAddressing =
        baseMesh().globalData().sharedPointAddr();

    // Global point data
    List< List< Map<vector> > > glData(Pstream::nProcs());
    forAll(glData, procI)
    {
        glData[procI] =
            List< Map<vector> >
            (
                baseMesh().globalData().nGlobalPoints(),
                Map<vector>()
            );
    }

    forAll(subMeshes(), meshI)
    {
        const labelList& pointMap = subMeshes()[meshI].pointMap();
        const vectorField& subMeshPointDI =
            subMeshFields[meshI].internalField();

        forAll(pointMap, pointI)
        {
            const label curMeshPoint = pointMap[pointI];
            const bool sharedPoint(findIndex(spLabels, curMeshPoint) != -1);

            if (sharedPoint)
            {
                const label k = findIndex(spLabels, curMeshPoint);
                const label curSpIndex = spAddressing[k];
                glData[Pstream::myProcNo()][curSpIndex].insert
                    (
                        meshI,
                        subMeshPointDI[pointI]
                    );
            }
            else
            {
                baseMeshFieldI[curMeshPoint] +=
                    subMeshPointDI[pointI]/noMat[curMeshPoint];
            }
        }
    }

    Pstream::gatherList(glData);
    Pstream::scatterList(glData);

    const int nGlobalPoints = baseMesh().globalData().nGlobalPoints();

    if (nGlobalPoints)
    {
        const int nSubMeshes = subMeshes().size();

        for (label k = 0; k < nGlobalPoints; k++)
        {
            // Current shared point index
            const label curSpIndex = findIndex(spAddressing, k);

            if (curSpIndex != -1)
            {
                List<label> matN(subMeshes().size(), 0);
                List<vector> matAvg(subMeshes().size(), vector::zero);

                forAll(glData, procI)
                {
                    const Map<vector>& curProcGlData = glData[procI][k];

                    for (label i = 0; i < nSubMeshes; i++)
                    {
                        if (curProcGlData.found(i))
                        {
                            matAvg[i] += curProcGlData[i];
                            matN[i]++;
                        }
                    }
                }

                label nMat = 0;
                vector avg = vector::zero;

                forAll(matAvg, matI)
                {
                    if (matN[matI])
                    {
                        matAvg[matI] /= matN[matI];
                        avg += matAvg[matI];
                        nMat++;
                    }
                }
                avg /= nMat;

                baseMeshFieldI[spLabels[curSpIndex]] = avg;
            }
        }
    }


    const int nIsolatedInterfacePoints =
        returnReduce(isolatedInterfacePoints().size(), sumOp<int>());

    if (nIsolatedInterfacePoints)
    {
        // Correct isolated points
        forAll(baseMeshField.boundaryField(), patchI)
        {
            if
            (
                isA<processorPointPatchVectorField>
                (
                    baseMeshField.boundaryField()[patchI]
                )
            )
            {
                const processorPointPatchVectorField& procPatchDispl =
                    dynamicCast<const processorPointPatchVectorField>
                    (
                        baseMeshField.boundaryField()[patchI]
                    );

                const processorPolyPatch& procPatch =
                    dynamicCast<const processorPolyPatch>
                    (
                        baseMesh().boundaryMesh()[patchI]
                    );

                vectorField pif
                (
                    procPatchDispl.patchInternalField
                    (
                        baseMeshField.internalField()
                    )
                );

                OPstream::write
                (
                    Pstream::commsTypes::blocking,
                    procPatch.neighbProcNo(),
                    reinterpret_cast<const char*>(pif.begin()),
                    pif.byteSize()
                );
            }
        }

        forAll(baseMeshField.boundaryField(), patchI)
        {
            if
            (
                baseMeshField.boundaryField()[patchI].type()
             == processorPointPatchVectorField::typeName
            )
            {
                if (Pstream::parRun())
                {
                    const processorPointPatchVectorField& procPatchDispl =
                        dynamicCast<const processorPointPatchVectorField>
                        (
                            baseMeshField.boundaryField()[patchI]
                        );

                    const processorPolyPatch& procPatch =
                        dynamicCast<const processorPolyPatch>
                        (
                            baseMesh().boundaryMesh()[patchI]
                        );

                    tmp<vectorField> tNgbProcPatchDispl
                    (
                        new vectorField(procPatchDispl.size(), vector::zero)
                    );
                    vectorField& ngbProcPatchDispl =
                        tNgbProcPatchDispl.ref();

                    IPstream::read
                    (
                        Pstream::commsTypes::blocking,
                        procPatch.neighbProcNo(),
                        reinterpret_cast<char*>(ngbProcPatchDispl.begin()),
                        ngbProcPatchDispl.byteSize()
                    );

                    const labelList& isoInterfacePoints =
                        isolatedInterfacePoints();
                    const labelListList& pointFaces = baseMesh().pointFaces();

                    forAll(isoInterfacePoints, pI)
                    {
                        const label curPoint = isoInterfacePoints[pI];

                        const labelList& curPointFaces = pointFaces[curPoint];

                        forAll(curPointFaces, fI)
                        {
                            const label faceID = curPointFaces[fI];
                            const label patchID =
                                baseMesh().boundaryMesh().whichPatch(faceID);

                            if (patchID == patchI)
                            {
                                const label curPatchPoint =
                                    baseMesh().boundaryMesh()
                                    [
                                        patchI
                                    ].meshPointMap()[curPoint];

                                baseMeshFieldI[curPoint] =
                                    ngbProcPatchDispl[curPatchPoint];

                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    baseMeshField.correctBoundaryConditions();

    // Re-evaluate componentMixed patches
//     forAll(baseMeshField.boundaryField(), patchI)
//     {
//         if
//         (
//             baseMeshField.boundaryField()[patchI].type()
//          == componentMixedPointPatchVectorField::typeName
//         )
//         {
//             baseMeshField.boundaryField()[patchI].evaluate();
//         }
//     }
}


// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //


template<class Type>
Foam::tmp< Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::solidSubMeshes::lookupBaseMeshVolField
(
    const word& fieldName,
    const fvMesh& subMesh
) const
{
    // Lookup the field from the base mesh
    const GeometricField<Type, fvPatchField, volMesh>& baseField =
        baseMesh().lookupObject< GeometricField<Type, fvPatchField, volMesh> >
        (
            fieldName
        );

    if (baseMesh().name() == subMesh.name())
    {
        // If the subMesh is the baseMesh then return the baseField
        return
            tmp< GeometricField<Type, fvPatchField, volMesh> >
            (
                new GeometricField<Type, fvPatchField, volMesh>
                (
                    fieldName + "_copy", baseField
                )
            );
    }
    else
    {
        // Find the fvMeshSubset corresponding to subMesh
        label curSubMeshID = -1;
        forAll(subMeshes(), meshI)
        {
            if (subMeshes()[meshI].subMesh().name() == subMesh.name())
            {
                curSubMeshID = meshI;
                break;
            }
        }

        if (curSubMeshID == -1)
        {
            FatalErrorIn
                (
                    "template<class Type>\n"
                    "Foam::tmp< Foam::GeometricField"
                    "<Type, Foam::fvPatchField, Foam::volMesh> >\n"
                    "Foam::solidSubMeshes::lookupBaseMeshVolField\n"
                    "(\n"
                    "    const word& fieldName,\n"
                    "    const fvMesh& subMesh\n"
                    ") const"
                )   << "SubMesh not found when looking for a field in the base "
                    << "mesh" << abort(FatalError);
        }

        // Return the baseField interpolated to the subMesh
        return subMeshes()[curSubMeshID].interpolate(baseField);
    }
}


template<class Type>
Foam::tmp< Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::solidSubMeshes::lookupBaseMeshSurfaceField
(
    const word& fieldName,
    const fvMesh& subMesh
) const
{
    // Lookup the field from the base mesh
    const GeometricField<Type, fvsPatchField, surfaceMesh>& baseField =
        baseMesh().lookupObject
        <
            GeometricField<Type, fvsPatchField, surfaceMesh>
        >
        (
            fieldName
        );

    if (baseMesh().name() == subMesh.name())
    {
        // If the subMesh is the baseMesh then return the baseField
        return
            tmp< GeometricField<Type, fvsPatchField, surfaceMesh> >
            (
                new GeometricField<Type, fvsPatchField, surfaceMesh>
                (
                    fieldName + "_copy", baseField
                )
            );
    }
    else
    {
        // Find the fvMeshSubset corresponding to subMesh
        label curSubMeshID = -1;
        forAll(subMeshes(), meshI)
        {
            if (subMeshes()[meshI].subMesh().name() == subMesh.name())
            {
                curSubMeshID = meshI;
                break;
            }
        }

        if (curSubMeshID == -1)
        {
            FatalErrorIn
                (
                    "template<class Type>\n"
                    "Foam::tmp< Foam::GeometricField"
                    "<Type, Foam::fvsPatchField, Foam::surfaceMesh> >\n"
                    "Foam::solidSubMeshes::lookupBaseMeshSurfaceField\n"
                    "(\n"
                    "    const word& fieldName,\n"
                    "    const fvMesh& subMesh\n"
                    ") const"
                )   << "SubMesh not found when looking for a field in the base "
                    << "mesh" << abort(FatalError);
        }

        // Return the baseField interpolated to the subMesh
        return subMeshes()[curSubMeshID].interpolate(baseField);
    }
}


// ************************************************************************* //
