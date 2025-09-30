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

Class
    Foam::cellRemovalFvMesh

Description
    Mesh allowing cell removal based on a run-time selectable cell removal law.

SourceFiles
    cellRemovalFvMesh.C

Author
    Philip Cardiff, UCD/UT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef cellRemovalFvMesh_H
#define cellRemovalFvMesh_H

#include "dynamicFvMesh.H"
#include "surfaceFields.H"
#include "mapPolyMesh.H"
#include "cellRemovalLaw.H"

#include "fvMeshSubset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cellRemovalFvMesh Declaration
\*---------------------------------------------------------------------------*/

class cellRemovalFvMesh
:
    public dynamicFvMesh
{
protected:
    // Protected Data Members

        //- Settings dictionary
        dictionary dict_;

        //- Remove dead cells: these are small islands of cells which have
        //  become disconnected from the main mesh. We define an island as a
        //  bunch of cells disconnected from the main mesh having less than 1%
        //  of the total number of cells in the mesh
        const Switch removeDeadCells_;

        //- Law for governing which cells are removed
        autoPtr<cellRemovalLaw> lawPtr_;

        //- Current time index
        label curIndex_;

        //- Flag to save the subset cells
        mutable bool saveSubset_;

        //- Removed subMesh
        mutable autoPtr<fvMeshSubset> subsetter_;


    // Private Member Functions

        //- Compute the map from face index to patch index for exposed faces
        labelList patchFaceMap
        (
            const mapPolyMesh& mpm,
            const labelList& exposedFaces
        ) const;

        //- Mark dead cells
        label addDeadCells(labelList& cellsToRemove);

        //- Sync faces across coupled patches
        void syncVisitedFaces
        (
            PackedBoolList& visited,
            PackedBoolList& marked,
            labelHashSet& newFacesToVisit
        ) const;

        //- Update volField values on the newly exposed faces
        //  The values on the new faces are extrapolated from the internal
        //  field; currently we use zero gradient extrapolation
        template<class Type>
        void updateVolFieldsExposedFaces
        (
            const mapPolyMesh& map,
            const labelList& exposedFaces,
            const labelList& patchFaceMap
        ) const;

        //- Save geometric fields
        template<class Type, template<class> class Patch, class Mesh>
        void saveGeoFields(objectRegistry&) const;

        //- Disallow default bitwise copy construct
        cellRemovalFvMesh(const cellRemovalFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const cellRemovalFvMesh&);


public:

    //- Runtime type information
    TypeName("cellRemovalFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit cellRemovalFvMesh(const IOobject& io);


    // Destructor

        virtual ~cellRemovalFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Apply topological changes
        virtual bool refine(const bool correctError = true)
        {
            return false;
        }

        //- Mark to save the previous regions when removing cells
        void subsetRemovedCells() const
        {
            saveSubset_ = true;
        }

        //- Return if cells have been removed
        bool active() const
        {
            return subsetter_.valid();
        }

        //- Return the subset region
        const fvMesh& removedMesh() const
        {
            return subsetter_->subMesh();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "cellRemovalFvMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
