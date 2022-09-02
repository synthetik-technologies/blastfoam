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

#include "cellRemovalFvMesh.H"
#include "fvPatchField.H"
#include "calculatedFvPatchField.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::cellRemovalFvMesh::updateVolFieldsExposedFaces
(
    const mapPolyMesh& map,
    const labelList& exposedFaces
) const
{
    // Get reverse face map
    const labelList& revFaceMap = map.reverseFaceMap();

    // Create a hash set from the expoed facess for easy searching
    //const labelHashSet exposedFacesSet(exposedFaces);

    // Read volField objects from object registry
    HashTable<const GeometricField<Type, fvPatchField, volMesh>*> fields
    (
        thisDb().objectRegistry::template lookupClass
        <GeometricField<Type, fvPatchField, volMesh> >()
    );

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

        Field<Type>& fieldI = field.primitiveFieldRef();

        if (debug)
        {
            Info<< "    volField " << fieldIter()->name() << endl;
        }

        // Check if there is a surface field by the same name suffixed with 'f'
        bool surfaceFieldExists = false;
        const GeometricField<Type, fvsPatchField, surfaceMesh>*
            surfaceFieldPtr = NULL;
        if
        (
            foundObject<GeometricField<Type, fvsPatchField, surfaceMesh> >
            (
                "interpolate(" + field.name() + ')'
            )
        )
        {
            surfaceFieldExists = true;

            surfaceFieldPtr =
              &(
                    lookupObject
                    <
                        GeometricField<Type, fvsPatchField, surfaceMesh>
                    >("interpolate(" + field.name() + ')')
               );

            if (debug)
            {
                Info<< "    surfaceField " << surfaceFieldPtr->name() << endl;
            }
        }

        // Initialise field on newly exposed faces

        forAll(exposedFaces, fI)
        {
            // Get new face ID
            label newFaceID = revFaceMap[exposedFaces[fI]];

            // Find the patch ID
            const label patchID = boundaryMesh().whichPatch(newFaceID);

            if (patchID == -1)
            {
                FatalErrorIn
                (
                    "void Foam::cellRemovalFvMesh::updateVolFieldsExposedFaces"
                )   << "exposed face is not on the boundary!? What's going on?"
                    << abort(FatalError);
            }

            const label start = boundaryMesh()[patchID].start();

            // Get local face ID
            const label newLocalFaceID = newFaceID - start;

            // Get face cell ID
            const label faceCellID =
                boundaryMesh()[patchID].faceCells()[newLocalFaceID];

            // Set the new face value to be the face cell value
            field.boundaryFieldRef()[patchID][newLocalFaceID] =
                fieldI[faceCellID];

            if (surfaceFieldExists)
            {
                GeometricField<Type, fvsPatchField, surfaceMesh>&
                    surfaceField =
                    const_cast
                    <
                        GeometricField<Type, fvsPatchField, surfaceMesh>&
                    >(*surfaceFieldPtr);

                // Set the new surface face value to be the face cell value
                surfaceField.boundaryFieldRef()[patchID][newLocalFaceID] =
                    fieldI[faceCellID];
            }
        }
    }
}


// ************************************************************************* //
