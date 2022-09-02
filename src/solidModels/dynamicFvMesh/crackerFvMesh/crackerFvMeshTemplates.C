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

#include "crackerFvMesh.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::crackerFvMesh::globalCrackField
(
    const Field<Type>& localCrackField
) const
{
    tmp<Field<Type> > tGlobalCrackField
    (
        new Field<Type>(this->globalCrackSize(), pTraits<Type>::zero)
    );

    label globalIndex = this->localCrackStart();

    label localCrackSize = boundaryMesh()[crackPatchID_.index()].size();

    for (label i=0; i<localCrackSize; i++)
    {
        tGlobalCrackField.ref()[globalIndex++] = localCrackField[i];
    }

    reduce(tGlobalCrackField.ref(), sumOp<List<Type> >());

    return tGlobalCrackField;
}


template<class Type>
void Foam::crackerFvMesh::updateVolFieldBrokenFaces
(
    const labelList& faceMap,
    const labelList& facesToBreak,
    const labelList& coupledFacesToBreak
) const
{
    // Create a hash set from the list of faces for easy searching
    const labelHashSet facesToBreakSet(facesToBreak);
    const labelHashSet coupledFacesToBreakSet(coupledFacesToBreak);

    // Cast the mesh to a crackerFvMesh

    const crackerFvMesh& mesh = refCast<const crackerFvMesh>(*this);

    // Read volField objects from object registry
    HashTable<const GeometricField<Type, fvPatchField, volMesh>*> fields
    (
        mesh.thisDb().objectRegistry::template lookupClass
       <GeometricField<Type, fvPatchField, volMesh> >()
    );


    const label cohesivePatchID = crackPatchID_.index();
    const label start = mesh.boundaryMesh()[cohesivePatchID].start();
    const label cohesivePatchSize = mesh.boundaryMesh()[cohesivePatchID].size();

    for
    (
        typename HashTable<const GeometricField<Type, fvPatchField, volMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        // Read field
        GeometricField<Type, fvPatchField, volMesh>& field =
            const_cast<GeometricField<Type, fvPatchField, volMesh>&>
            (*fieldIter());

        if (debug)
        {
            Info<< "    volField " << fieldIter()->name() << endl;
        }

        // Local crack field
        Field<Type> fieldpI
        (
            field.boundaryField()[cohesivePatchID].patchInternalField()
        );

        // Global crack fields
        Field<Type> gFieldpI(mesh.globalCrackField(fieldpI));

        // Check if there is a surface equivalent to the vol field
        bool surfaceFieldExists = false;
        const GeometricField<Type, fvsPatchField, surfaceMesh>*
            surfaceFieldPtr = NULL;
        if
        (
            mesh.foundObject<GeometricField<Type, fvsPatchField, surfaceMesh> >
            (
                "interpolate(" + field.name() + ')'
            )
        )
        {
            surfaceFieldExists = true;

            surfaceFieldPtr =
              &(
                    mesh.lookupObject
                    <
                        GeometricField<Type, fvsPatchField, surfaceMesh>
                    >("interpolate(" + field.name() + ')')
               );

            if (debug)
            {
                Info<< "    surfaceField " << surfaceFieldPtr->name() << endl;
            }
        }

        // Initialise field on new cohesive face

        const labelList& gcfa = mesh.globalCrackFaceAddressing();

        label globalIndex = mesh.localCrackStart();

        for (label i = 0; i < cohesivePatchSize; i++)
        {
            label oldFaceIndex = faceMap[start + i];

            if
            (
                facesToBreakSet.found(oldFaceIndex)
             || coupledFacesToBreakSet.found(oldFaceIndex)
            )
            {
                // Set the new face value to be the average of the previously
                // adjoining cells

                field.boundaryFieldRef()[cohesivePatchID][i] =
                    0.5
                   *(
                        gFieldpI[globalIndex] + gFieldpI[gcfa[globalIndex]]
                    );

                if (surfaceFieldExists)
                {
                    GeometricField<Type, fvsPatchField, surfaceMesh>&
                        surfaceField =
                        const_cast
                        <
                            GeometricField<Type, fvsPatchField, surfaceMesh>&
                        >(*surfaceFieldPtr);

                    surfaceField.boundaryFieldRef()[cohesivePatchID][i] =
                      field.boundaryField()[cohesivePatchID][i];
                }

                globalIndex++;
            }
            else
            {
                globalIndex++;
            }
        }
    }
}


// ************************************************************************* //
