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
    const labelList& exposedFaces,
    const labelList& patchFaceMap
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolTypeField;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfTypeField;

    // Get reverse face map
    const labelList& revFaceMap = map.reverseFaceMap();

    // Create a hash set from the expoed facess for easy searching
    //const labelHashSet exposedFacesSet(exposedFaces);

    // Read volField objects from object registry
    HashTable<const VolTypeField*> fields
    (
        thisDb().template lookupClass<VolTypeField>()
    );

    forAllConstIter
    (
        typename HashTable<const VolTypeField*>,
        fields,
        fieldIter
    )
    {
        DebugInfo<< "Updating volField exposed faces for " << fieldIter()->name() << endl;

        // Read field
        VolTypeField& field = const_cast<VolTypeField&>(*fieldIter());

        // Check if there is a surface field by the same name suffixed with 'f'
        const word ffieldName("interpolate(" + field.name() + ")");
        if (foundObject<SurfTypeField>(ffieldName))
        {

            SurfTypeField& surfaceField =
                lookupObjectRef<SurfTypeField>(ffieldName);

            DebugInfo<< "    Using " << surfaceField.name()
                << " for correction" << endl;

            typename VolTypeField::Boundary& bfield = field.boundaryFieldRef();
            typename SurfTypeField::Boundary& bsurfaceField =
                surfaceField.boundaryFieldRef();

            // Initialise field on newly exposed faces
            forAll(exposedFaces, fi)
            {
                const label oldFaceID = exposedFaces[fi];
                // Get new face ID
                label newFaceID = revFaceMap[oldFaceID];

                // Find the patch ID
                const label patchID = patchFaceMap[fi];

                if (patchID == -1)
                {
                    FatalErrorInFunction
                        << "exposed face is not on the boundary!? What's going on?"
                        << abort(FatalError);
                }

                const label start = boundaryMesh()[patchID].start();

                // Get local face ID
                const label newLocalFaceID = newFaceID - start;

                // Get face cell ID
                const label faceCellID =
                    boundaryMesh()[patchID].faceCells()[newLocalFaceID];

                // Set the new face value to be the previous face value
                bfield[patchID][newLocalFaceID] = field[faceCellID];
                bsurfaceField[patchID][newLocalFaceID] = field[faceCellID];
            }
        }
        else
        {
            typename VolTypeField::Boundary& bfield = field.boundaryFieldRef();

            // Initialise field on newly exposed faces
            forAll(exposedFaces, fi)
            {
                // Get new face ID
                label newFaceID = revFaceMap[exposedFaces[fi]];

                // Find the patch ID
                const label patchID = patchFaceMap[fi];

                const label start = boundaryMesh()[patchID].start();

                // Get local face ID
                const label newLocalFaceID = newFaceID - start;

                // Get face cell ID
                const label faceCellID =
                    boundaryMesh()[patchID].faceCells()[newLocalFaceID];

                // Set the new face value to be the face cell value
                bfield[patchID][newLocalFaceID] = field[faceCellID];
            }
        }
    }
}


template<class Type, template<class> class Patch, class Mesh>
void Foam::cellRemovalFvMesh::saveGeoFields
(
    objectRegistry& obr
) const
{
    typedef GeometricField<Type, Patch, Mesh> GeoField;
    HashTable<GeoField*> fields = obr.lookupClass<GeoField>();

    forAllIter
    (
        typename HashTable<GeoField*>,
        fields,
        iter
    )
    {
        const GeoField& field = *iter();
        GeoField* interpField = subsetter_->interpolate(*iter()).ptr();
        interpField->rename(field.name());
        interpField->writeOpt() = field.writeOpt();
        interpField->store(interpField);
    }
}


// ************************************************************************* //
