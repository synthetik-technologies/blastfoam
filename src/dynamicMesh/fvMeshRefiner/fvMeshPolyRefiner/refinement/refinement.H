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

Class
    Foam::refinement

Description
    Abstract base class for adaptive mesh refinement using the mesh modifier
    engine. The class provides common interface and functionalities for 3D
    polyhedral refinement and 2D prismatic refinement.

    The common interface includes (pure virtuals) following member functions:
    - setCellsToRefine
    - setSplitPointsToUnrefine

    With a lot of ordinary protected member functions which are used by both
    derived classes.

    Note: I've written it this way in order to avoid unnecesasry code
    duplication, but I'm 99% sure that if someone else wants to write additional
    refinement strategy (e.g. directional refinement) as derived class, the
    interface will need to change.

SourceFiles
    refinement.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

Notes
    Generalisation of hexRef8 for polyhedral cells and refactorisation using
    polyMesh modifier engine.

\*---------------------------------------------------------------------------*/

#ifndef refinement_H
#define refinement_H

#include "Time.H"
#include "polyMesh.H"
#include "labelIOField.H"
#include "removeFaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class refinement Declaration
\*---------------------------------------------------------------------------*/

class refinement
{
protected:

    // Protected data

        //- Reference to polyMesh for easy access in helper functions
        const polyMesh& mesh_;


        // Refinement control and handling

            //- Cell refinement level
            mutable labelIOList cellLevel_;

            //- Point refinement level
            mutable labelIOList pointLevel_;

            //- Face remover engine
            mutable removeFaces faceRemover_;

            //- Switch whether to use edge based consistency on refinement
            bool edgeBasedConsistency_;


    // Protected Pure Virtual Member Functions

        // Global topology modification functions (operate on whole polyMesh)

            //- Set refinement instruction
            virtual void setRefinement
            (
                polyTopoChange& meshMod,
                const labelList& cellsToRefine
            ) const = 0;

            //- Set unrefinement instruction
            virtual void setUnrefinement
            (
                polyTopoChange& meshMod,
                const labelList& splitPointsToUnrefine
            ) const = 0;


    // Protected Member Functions

        // Useful helper functions used by derived classes

            //- Set file instance for cellLevel_ and pointLevel_
            void setInstance(const fileName& inst) const;


        // Local topology modification functions (operate on cells/faces)

            //- Adds a face on top of existing faceI. Reverses if nessecary
            label addFace
            (
                polyTopoChange& ref,
                const label faceI,
                const face& newFace,
                const label own,
                const label nei
            ) const;

            //- Adds internal face from point. No checks on reversal
            label addInternalFace
            (
                polyTopoChange& ref,
                const label meshFaceI,
                const label meshPointI,
                const face& newFace,
                const label own,
                const label nei
            ) const;

            //- Modifies existing faceI for either new owner/neighbour or new
            //  face points. Reverses if nessecary
            void modifyFace
            (
                polyTopoChange& ref,
                const label faceI,
                const face& newFace,
                const label own,
                const label nei
            ) const;


        // Topological change helper functions

            //- Store vertices from startFp up to face split point.
            //  Used when splitting face into n faces where n is the number of
            //  points in a face (or number of edges)
            void walkFaceToMid
            (
                const labelList& edgeMidPoint,
                const label cLevel,
                const label faceI,
                const label startFp,
                DynamicList<label>& faceVerts
            ) const;

            //- Same as walkFaceToMid but now walk back
            void walkFaceFromMid
            (
                const labelList& edgeMidPoint,
                const label cLevel,
                const label faceI,
                const label startFp,
                DynamicList<label>& faceVerts
            ) const;

            //- Get index of point with minimum point level
            label findMinLevel(const labelList& f) const;

            //- Get index of point with maximum point level
            label findMaxLevel(const labelList& f) const;

            //- Count number of vertices <= anchorLevel for a given face
            label countAnchors
            (
                const labelList& f,
                const label anchorLevel
            ) const;


        // Debug functions

            //- Check orientation of added internal face
            void checkInternalOrientation
            (
                polyTopoChange& ref,
                const label cellI,
                const label faceI,
                const point& ownPt,
                const point& neiPt,
                const face& newFace
            ) const;

            //- Check orientation of a new boundary face
            void checkBoundaryOrientation
            (
                polyTopoChange& ref,
                const label cellI,
                const label faceI,
                const point& ownPt,
                const point& boundaryPt,
                const face& newFace
            ) const;


        // Refinement/unrefinement consistency checks

            //- Updates cellsToRefine such that a face consistent 2:1 refinement
            //  is obtained. Returns local number of cells changed
            label faceConsistentRefinement
            (
                PackedBoolList& refineCell,
                const bool maxSet
            ) const;

            //- Updates cellsToRefine such that an edge consistent 4:1 refinement
            //  is obtained. Returns local number of cells changed
            label edgeConsistentRefinement
            (
                PackedBoolList& refineCell,
                const bool maxSet
            ) const;

            //- Updates cellsToUnrefine such that a face consistent 2:1
            //  unrefinement is obtained. Returns local number of cells changed
            label faceConsistentUnrefinement
            (
                PackedBoolList& refineCell,
                const bool maxSet
            ) const;

            //- Updates cellsToUnrefine such that an edge consistent 4:1
            //  unrefinement is obtained. Returns local number of cells changed
            label edgeConsistentUnrefinement
            (
                PackedBoolList& refineCell,
                const bool maxSet
            ) const;


        // Copy control

            //- Disallow default bitwise copy construct
            refinement(const refinement&);

            //- Disallow default bitwise assignment
            void operator=(const refinement&);


public:

    //- Runtime type information
    TypeName("refinement");


    // Constructors

        //- Construct from dictionary
        refinement
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const bool read = true
        );


    //- Destructor
    virtual ~refinement();


    // Member Functions

        // Access

            //- Return refinement cell level
            const labelIOList& cellLevel() const
            {
                return cellLevel_;
            }

            //- Return refinement point level
            const labelIOList& pointLevel() const
            {
                return pointLevel_;
            }


        // Edit

            //- Set cells to refine given a list of refinement
            //  candidates. Refinement candidates are extended within the
            //  function due to possible 4:1 conflicts and specified number of
            //  buffer layers.
            virtual labelList consistentRefinement
            (
                const labelList& refinementCellCandidates,
                const bool maxSet
            ) const;

            //- Set split points to unrefine given a list of all mesh points
            //  that are candidates for unrefinement. Split points are
            //  determined as a subset of unrefinement candidates, avoiding
            //  splitting points of cells that are going to be refined at the
            //  same time and ensuring consistent unrefinement.
            //  Note: must be called AFTER setCellsToRefine
            virtual labelList consistentUnrefinement
            (
                const labelList& unrefinementPointCandidates,
                const bool maxSet
            ) const = 0;


        // Inherited interface from polyMeshModifier


            //- Insert the polyhedral refinement into the
            //  topological change and return the cells that were refined
            //  (new indices)
            virtual autoPtr<mapPolyMesh> refine
            (
                polyMesh& mesh,
                const labelList& cellsToRefine
            );

            //- Insert the polyhedral refinement/unrefinement into the
            //  topological change
            virtual bool unrefine
            (
                polyMesh& mesh,
                const labelList& unrefinementPointCandidates
            );

            //- Force recalculation of locally stored data on topological change
            virtual void updateMesh(const mapPolyMesh&);

            //- Distribute additional information
            virtual void distribute(const mapDistributePolyMesh& map);


        //- Write
        virtual bool write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
