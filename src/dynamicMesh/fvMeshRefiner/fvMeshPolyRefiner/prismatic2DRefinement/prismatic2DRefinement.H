/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    Foam::prismatic2DRefinement

Description
    Isotropic refinement of prismatic cells in 2D using the mesh modifier
    engine. Used for 2D cases instead of polyhedralRefinement which carries
    unnecessary overhead in terms of number of cells since it splits the cell in
    all directions.

    Each prismatic cell is split by the following procedure:
    1. Adding points at the face centres and edge centres of all faces found on
       an special patch: empty or wedge.
    2. Adding n cells per existing cell where n is the number of corner points
       at the face on special patch (empty or wedge).
    3. Splitting each of the faces on special patch (empty or wedge) into
       multiple faces going from: existing corner point -> new edge centre point
       -> new face centre point -> other new edge centre point (sharing the same
       corner point)
    4. Spliiting each of the faces not on an special patch (empty or wedge) into
       two faces going from: existing corner point -> existing corner point on
       the other side -> new edge centre point on the other side -> new edge
       centre point on my side
    4. Adding internal faces going from:
       new edge centre point -> new face centre point -> new other face
       centre point on the other side -> new other edge mid point on the other
       side

    It is an error to try and run this on anything except a 2D mesh.

SourceFiles
    prismatic2DRefinement.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef prismatic2DRefinement_H
#define prismatic2DRefinement_H

#include "refinement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class prismatic2DRefinement Declaration
\*---------------------------------------------------------------------------*/

class prismatic2DRefinement
:
    public refinement
{
private:

    // Private Member Functions

        // Helper functions

            //- Get least cell level such that the face has at least nPoints
            //  points smaller than the level
            label getAnchorLevel
            (
                const label faceI,
                const label nPoints
            ) const;

            //- Append given face into a dynamic list containing split faces
            //  that will be split into two faces (third parameter). Additionaly
            //  append information on which of the two edges of the face are on
            //  special patch into a dynamic list (fourth parameter)
            void appendFaceSplitInfo
            (
                const label& faceI,
                const boolList& edgeOnEmptyPatch,
                const labelList& edgeMidPoint,
                DynamicList<label>& splitFacesIntoTwo,
                DynamicList<Pair<label> >& splitFacesEmptyEdges
            ) const;


        // Topological change helper functions

            //- Set new owner and neighbour given anchor pointI, faceI and the
            //  necessary mapping
            void setNewFaceNeighbours
            (
                const HashTable
                <
                    label,
                    Pair<label>,
                    Hash<FixedList<label, 2> >
                >& pointCellToAddedCellMap,
                const labelListList& cellAddedCells,
                const label& faceI,
                const label& pointI,

                label& own,
                label& nei
            ) const;

            //- Get index of point with minimum point level of a face across two
            //  connected edges starting from a local point index.
            //  Example: starting from point with level 1 in the upper left
            //  corner, finds point index of the point with level 0 which is on
            //  the same face, connected with edge to original point
            //  1------1
            //  |
            //  |
            //  0
            //  Note: passing face edges and mesh edges as parameters to avoid
            //  fetching them from mesh due to lazy evaluation
            label findMinEdgeConnectedLevel
            (
                const label& fpI,
                const label& faceI,
                const face& f,
                const labelList& fEdges,
                const edgeList& meshEdges
            ) const;

            //- Store two face mids when adding internal faces
            void addFaceMids
            (
                const labelList& faceMidPoint,
                const boolList& faceOnEmptyPatch,
                const label& faceI,
                const label& cellI,
                face& newFace
            ) const;


        // Debug functions

            //- Check orientation of a split face
            void checkNewFaceOrientation
            (
                polyTopoChange& ref,
                const label& faceI,
                const face& newFace
            ) const;


        // Copy control

            //- Disallow default bitwise copy construct
            prismatic2DRefinement(const prismatic2DRefinement&);

            //- Disallow default bitwise assignment
            void operator=(const prismatic2DRefinement&);


protected:

    // Protected Pure Virtual Member Functions

        // Global topology modification functions (operate on whole polyMesh)

            //- Set refinement instruction
            virtual void setRefinement
            (
                polyTopoChange& ref,
                const labelList& cellsToRefine
            ) const;

            //- Set unrefinement instruction
            virtual void setUnrefinement
            (
                polyTopoChange& ref,
                const labelList& splitPointsToUnrefine
            ) const;


public:

    //- Runtime type information
    TypeName("prismatic2DRefinement");


    // Constructors

        //- Construct from dictionary
        prismatic2DRefinement
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const bool read = true
        );


    //- Destructor
    virtual ~prismatic2DRefinement();


    // Member Functions

        // Edit

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
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
