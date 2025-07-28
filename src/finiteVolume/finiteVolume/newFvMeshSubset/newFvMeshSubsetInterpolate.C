/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "newFvMeshSubset.H"
#include "calculatedFvPatchFields.H"
#include "calculatedFvsPatchFields.H"
#include "emptyFvPatchFields.H"
#include "emptyFvsPatchFields.H"
#include "calculatedPointPatchFields.H"

#include "symmetryFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > newFvMeshSubset::meshToMesh
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const fvMesh& sMesh,
    const labelList& patchMap,
    const labelList& cellMap,
    const labelList& faceMap
)
{
    // Create and map the internal-field values
    Field<Type> internalField(vf.internalField(), cellMap);

    // Create and map the patch field values
    PtrList<fvPatchField<Type> > patchFields(patchMap.size());

    forAll (patchFields, patchI)
    {
        // Set the first one by hand as it corresponds to the
        // exposed internal faces.  Additional interpolation can be put here
        // as necessary.  HJ, date deleted
        if (patchMap[patchI] == -1)
        {
            patchFields.set
            (
                patchI,
                new emptyFvPatchField<Type>
                (
                    sMesh.boundary()[patchI],
                    DimensionedField<Type, volMesh>::null()
                )
            );
        }
        else
        {
            // Construct addressing
            const fvPatch& subPatch = sMesh.boundary()[patchI];
            const fvPatch& basePatch = vf.mesh().boundary()[patchMap[patchI]];
            label baseStart = basePatch.patch().start();
            label baseSize = basePatch.size();

            labelList directAddressing(subPatch.size());

            forAll(directAddressing, i)
            {
                label baseFaceI = faceMap[subPatch.patch().start()+i];

                if (baseFaceI >= baseStart && baseFaceI < baseStart + baseSize)
                {
                    directAddressing[i] = baseFaceI-baseStart;
                }
                else
                {
                    // Mapped from internal face. Do what? Map from element
                    // 0 for now.
                    directAddressing[i] = 0;
                }
            }

            patchFields.set
            (
                patchI,
                fvPatchField<Type>::New
                (
                    vf.boundaryField()[patchMap[patchI]],
                    sMesh.boundary()[patchI],
                    DimensionedField<Type, volMesh>::null(),
                    patchFieldSubset
                    (
                        vf.boundaryField()[patchMap[patchI]].size(),
                        directAddressing
                    )
                )
            );

            // What to do with exposed internal faces if put into this patch?
        }
    }


    // Create the complete field from the pieces
    tmp<GeometricField<Type, fvPatchField, volMesh> > tresF
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "subset"+vf.name(),
                sMesh.time().timeName(),
                sMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            sMesh,
            vf.dimensions(),
            internalField,
            patchFields
        )
    );

    return tresF;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > newFvMeshSubset::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Get reference to the subset mesh
    const fvMesh& sMesh = subMesh();

    // Create and map the internal-field values
    Field<Type> internalField(vf.internalField(), cellMap());

    // Create and map the patch field values
    const labelList& pm = patchMap();

    // Create and map the patch field values
    PtrList<fvPatchField<Type> > patchFields(pm.size());

    label internalFacesPatchIndex = -1;

    forAll (patchFields, patchI)
    {
        // Set the first one by hand as it corresponds to the
        // exposed internal faces.  Additional interpolation can be put here
        // as necessary.  HJ, date deleted
        if (pm[patchI] == -1)
        {
            // Bug fix. Zeljko Tukovic, 10/Mar/2010
            internalFacesPatchIndex = patchI;

            patchFields.set
            (
                patchI,
                new calculatedFvPatchField<Type>
                (
                    sMesh.boundary()[patchI],
                    DimensionedField<Type, volMesh>::null()
                )
            );
        }
        else
        {
            patchFields.set
            (
                patchI,
                fvPatchField<Type>::New
                (
                    vf.boundaryField()[pm[patchI]],
                    sMesh.boundary()[patchI],
                    DimensionedField<Type, volMesh>::null(),
                    patchFieldSubset(*this, patchI)
                )
            );
        }
    }

    // Linear interpolation for last patch
    if (internalFacesPatchIndex > -1)
    {
        const Field<Type>& vfI = vf.internalField();
        const scalarField& w = baseMesh().weights().internalField();
        const labelList& own = baseMesh().faceOwner();
        const labelList& ngb = baseMesh().faceNeighbour();

        Field<Type>& lastPatchField = patchFields[internalFacesPatchIndex];

        label lastPatchStart =
            sMesh.boundaryMesh()[internalFacesPatchIndex].start();

        const labelList& fm = faceMap();

        forAll(lastPatchField, faceI)
        {
            if (fm[lastPatchStart + faceI] < baseMesh().nInternalFaces())
            {
                lastPatchField[faceI] =
                    w[fm[lastPatchStart + faceI]]*
                    vfI[own[fm[lastPatchStart + faceI]]]
                  + (1.0 - w[fm[lastPatchStart + faceI]])*
                    vfI[ngb[fm[lastPatchStart + faceI]]];
            }
            else
            {
                label patchID =
                    baseMesh().boundaryMesh().whichPatch
                    (
                        fm[lastPatchStart + faceI]
                    );

                label localFaceIndex =
                    fm[lastPatchStart + faceI]
                  - baseMesh().boundaryMesh()[patchID].start();

                lastPatchField[faceI] =
                    vf.boundaryField()[patchID][localFaceIndex];
            }
        }
    }

    // Create the complete field from the pieces
    tmp<GeometricField<Type, fvPatchField, volMesh> > tresF
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "subset"+vf.name(),
                sMesh.time().timeName(),
                sMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            sMesh,
            vf.dimensions(),
            internalField,
            patchFields
        )
    );

    return tresF;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
newFvMeshSubset::meshToMesh
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& vf,
    const fvMesh& sMesh,
    const labelList& patchMap,
    const labelList& faceMap
)
{
    // Create and map the internal-field values
    Field<Type> internalField
    (
        vf.internalField(),
        SubList<label>
        (
            faceMap,
            sMesh.nInternalFaces()
        )
    );

    // Create and map the patch field values
    PtrList<fvsPatchField<Type> > patchFields(patchMap.size());

    forAll (patchFields, patchI)
    {
        // Set the first one by hand as it corresponds to the
        // exposed internal faces.  Additional interpolation can be put here
        // as necessary.
        if (patchMap[patchI] == -1)
        {
            patchFields.set
            (
                patchI,
                new emptyFvsPatchField<Type>
                (
                    sMesh.boundary()[patchI],
                    DimensionedField<Type, surfaceMesh>::null()
                )
            );
        }
        else
        {
            // Construct addressing
            const fvPatch& subPatch = sMesh.boundary()[patchI];
            const fvPatch& basePatch = vf.mesh().boundary()[patchMap[patchI]];
            label baseStart = basePatch.patch().start();
            label baseSize = basePatch.size();

            labelList directAddressing(subPatch.size());

            forAll(directAddressing, i)
            {
                label baseFaceI = faceMap[subPatch.patch().start()+i];

                if (baseFaceI >= baseStart && baseFaceI < baseStart+baseSize)
                {
                    directAddressing[i] = baseFaceI-baseStart;
                }
                else
                {
                    // Mapped from internal face. Do what? Map from element
                    // 0 for now.
                    directAddressing[i] = 0;
                }
            }

            patchFields.set
            (
                patchI,
                fvsPatchField<Type>::New
                (
                    vf.boundaryField()[patchMap[patchI]],
                    sMesh.boundary()[patchI],
                    DimensionedField<Type, surfaceMesh>::null(),
                    patchFieldSubset
                    (
                        vf.boundaryField()[patchMap[patchI]].size(),
                        directAddressing
                    )
                )
            );
        }
    }


    // Map exposed internal faces. Note: Only nessecary if exposed faces added
    // into existing patch but since we don't know that at this point...
    forAll(patchFields, patchI)
    {
        fvsPatchField<Type>& pfld = patchFields[patchI];

        label meshFaceI = pfld.patch().patch().start();

        forAll(pfld, i)
        {
            label oldFaceI = faceMap[meshFaceI++];

            if (oldFaceI < vf.internalField().size())
            {
                pfld[i] = vf.internalField()[oldFaceI];
            }
        }
    }

    // Create the complete field from the pieces
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tresF
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "subset"+vf.name(),
                sMesh.time().timeName(),
                sMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            sMesh,
            vf.dimensions(),
            internalField,
            patchFields
        )
    );

    return tresF;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
newFvMeshSubset::interpolate
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& vf
) const
{
    // Get reference to the subset mesh
    const fvMesh& sMesh = subMesh();

    // Create and map the internal-field values
    Field<Type> internalField
    (
        vf.internalField(),
        SubList<label>
        (
            faceMap(),
            sMesh.nInternalFaces()
        )
    );

    // Create and map the patch field values
    const labelList& pm = patchMap();

    // Create and map the patch field values
    PtrList<fvsPatchField<Type> > patchFields(pm.size());

    forAll (patchFields, patchI)
    {
        // Set the first one by hand as it corresponds to the
        // exposed internal faces.  Additional interpolation can be put here
        // as necessary.  HJ, date deleted
        if (pm[patchI] == -1)
        {
            patchFields.set
            (
                patchI,
                //calculatedFvsPatchField<Type> // PC: bugfix
                new calculatedFvsPatchField<Type>
                (
                    sMesh.boundary()[patchI],
                    DimensionedField<Type, surfaceMesh>::null()
                )
            );
        }
        else
        {
            patchFields.set
            (
                patchI,
                fvsPatchField<Type>::New
                (
                    vf.boundaryField()[pm[patchI]],
                    sMesh.boundary()[patchI],
                    DimensionedField<Type, surfaceMesh>::null(),
                    patchFieldSubset(*this, patchI)
                )
            );
        }
    }


    const labelList& fm = faceMap();

    // Map exposed internal faces. Note: Only nessecary if exposed faces added
    // into existing patch but since we don't know that at this point...
    forAll(patchFields, patchI)
    {
        fvsPatchField<Type>& pfld = patchFields[patchI];

        label meshFaceI = pfld.patch().patch().start();

        forAll(pfld, i)
        {
            label oldFaceI = fm[meshFaceI++];

            if (oldFaceI < vf.internalField().size())
            {
                pfld[i] = vf.internalField()[oldFaceI];
            }
        }
    }

    // Create the complete field from the pieces
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tresF
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "subset"+vf.name(),
                sMesh.time().timeName(),
                sMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            sMesh,
            vf.dimensions(),
            internalField,
            patchFields
        )
    );

    return tresF;
}


template<class Type>
tmp<GeometricField<Type, pointPatchField, pointMesh> >
newFvMeshSubset::interpolate
(
    const GeometricField<Type, pointPatchField, pointMesh>& vf
) const
{
    // Get reference to the subset mesh
    const pointMesh& sMesh = subPointMesh();

    // Create and map the internal-field values
    Field<Type> internalField(vf.internalField(), pointMap());

    // Create and map the patch field values
    const labelList& pm = patchMap();

    // Create and map the patch field values
    PtrList<pointPatchField<Type> > patchFields(pm.size());

    forAll (patchFields, patchI)
    {
        // Set the first one by hand as it corresponds to the
        // exposed internal faces.  Additional interpolation can be put here
        // as necessary.  HJ, date deleted
        if (pm[patchI] == -1)
        {
            patchFields.set
            (
                patchI,
                new CalculatedPointPatchField
                <
                    pointPatchField,
                    pointMesh,
                    pointPatch,
                    DummyMatrix,
                    Type
                >
                (
                    sMesh.boundary()[patchI],
                    DimensionedField<Type, pointMesh>::null()
                )
            );
        }
        else
        {
            // Construct addressing
            const pointPatch& basePatch =
                vf.mesh().boundary()[pm[patchI]];

            const labelList& meshPoints = basePatch.meshPoints();

            // Make addressing from mesh to patch point
            Map<label> meshPointMap(2*meshPoints.size());
            forAll(meshPoints, localI)
            {
                meshPointMap.insert(meshPoints[localI], localI);
            }

            // Find which subpatch points originate from which patch point
            const pointPatch& subPatch = sMesh.boundary()[patchI];
            const labelList& subMeshPoints = subPatch.meshPoints();

            // If mapped from outside patch use point 0 for lack of better.
            labelList directAddressing(subPatch.size(), 0);

            const labelList& ptMap = pointMap();

            forAll(subMeshPoints, localI)
            {
                // Get mesh point on original mesh.
                label meshPointI = ptMap[subMeshPoints[localI]];

                Map<label>::const_iterator iter =
                    meshPointMap.find(meshPointI);

                if (iter != meshPointMap.end())
                {
                    directAddressing[localI] = iter();
                }
            }

            patchFields.set
            (
                patchI,
                pointPatchField<Type>::New
                (
                    vf.boundaryField()[pm[patchI]],
                    subPatch,
                    DimensionedField<Type, pointMesh>::null(),
                    pointPatchFieldSubset
                    (
                        directAddressing
                    )
                )
            );
        }
    }

#ifndef OPENFOAMESIORFOUNDATION
    // Add the global patch
    if (isType<globalPointPatch>(sMesh.boundary()[sMesh.boundary().size()-1]))
    {
        patchFields.resize(pm.size() + 1);
        patchFields.set
        (
            pm.size(),
            new GlobalPointPatchField
            <
                pointPatchField,
                pointMesh,
                pointPatch,
                globalPointPatch,
                DummyMatrix,
                Type
            >
            (
                sMesh.boundary().globalPatch(),
                DimensionedField<Type, pointMesh>::null()
            )
        );
    }
#endif

    // Create the complete field from the pieces
    tmp<GeometricField<Type, pointPatchField, pointMesh> > tresF
    (
        new GeometricField<Type, pointPatchField, pointMesh>
        (
            IOobject
            (
                "subset"+vf.name(),
                vf.time().timeName(),
                subMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            sMesh,
            vf.dimensions(),
            internalField,
            patchFields
        )
    );

    return tresF;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
