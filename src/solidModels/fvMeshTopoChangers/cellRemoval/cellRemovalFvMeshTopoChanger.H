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
    Foam::fvMeshTopoChangers::cellRemoval

Description
    Mesh allowing cell removal based on a run-time selectable cell removal law.

SourceFiles
    cellRemovalFvMesh.C

Author
    Philip Cardiff, UCD/UT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef cellRemovalFvMeshTopoChanger_H
#define cellRemovalFvMeshTopoChanger_H

#include "fvMeshTopoChanger.H"
#include "surfaceFields.H"
#include "cellRemovalLaw.H"
#include "polyTopoChangeMap.H"
#include "fvMeshSubset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{

/*---------------------------------------------------------------------------*\
                           Class cellRemoval Declaration
\*---------------------------------------------------------------------------*/

class cellRemoval
:
    public fvMeshTopoChanger
{
protected:
    // Protected Data Members

        //- Remove dead cells: these are small islands of cells which have
        //  become disconnected from the main mesh. We define an island as a
        //  bunch of cells disconnected from the main mesh having less than 1%
        //  of the total number of cells in the mesh
        const Switch removeDeadCells_;

        //- Law for governing which cells are removed
        autoPtr<cellRemovalLaw> lawPtr_;

        //- Flag to save the subset cells
        mutable bool saveSubset_;

        //- Removed subMesh
        mutable autoPtr<fvMeshSubset> subsetter_;


    // Private Member Functions

        //- Compute the map from face index to patch index for exposed faces
        labelList patchFaceMap
        (
            const polyTopoChangeMap& mpm,
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
            const polyTopoChangeMap& map,
            const labelList& exposedFaces,
            const labelList& patchFaceMap
        ) const;

        //- Save geometric fields
        template<class Type, template<class> class Patch, class Mesh>
        void saveGeoFields(objectRegistry&) const;

        //- Disallow default bitwise copy construct
        cellRemoval(const cellRemoval&);

        //- Disallow default bitwise assignment
        void operator=(const cellRemoval&);


public:

    //- Runtime type information
    TypeName("cellRemoval");


    // Constructors

        //- Construct from IOobject
        explicit cellRemoval(fvMesh& mesh, const dictionary& dict);


    // Destructor

        virtual ~cellRemoval();


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

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshTopoChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "cellRemovalFvMeshTopoChangerTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
