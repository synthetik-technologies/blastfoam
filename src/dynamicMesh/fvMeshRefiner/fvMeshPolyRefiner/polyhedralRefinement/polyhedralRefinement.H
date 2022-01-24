/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
24-01-2022 Synthetik Applied    :   Added support for axis-symmetric
           Technologies                 cases and load balancing
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
    Foam::polyhedralRefinement

Description
    Isotropic refinement of polyhedral cells using the mesh modifier engine

    Each polyhedral cell is split by the following procedure:
    1. Adding points at the edge centre, face centre and cell centre,
    2. Adding cells n cells where n is the number of points of the cell,
    3. Splitting each face into multiple faces going from:
       existing corner point -> new edge centre point -> new face centre
       point -> other new edge centre point (sharing the same corner point)
    4. Adding internal faces going from:
       new edge centre point -> new face centre point -> new cell centre
       point -> other new face centre point (sharing the same edge)

SourceFiles
    polyhedralRefinement.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef polyhedralRefinement_H
#define polyhedralRefinement_H

#include "refinement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class polyhedralRefinement Declaration
\*---------------------------------------------------------------------------*/

class polyhedralRefinement
:
    public refinement
{
private:

    // Private Member Functions

        // Helper functions

            //- Get least cell level such that the face has at least three
            //  points smaller than the level
            label getAnchorLevel(const label faceI) const;


        // Local topology modification functions (operate on cells/faces)

            //- Create all internal faces of split cellI into n cells where n is the
            //  number of cell points
            void createInternalFaces
            (
                const labelListList& cellAnchorPoints,
                const labelListList& cellAddedCells,
                const labelList& cellMidPoint,
                const labelList& faceMidPoint,
                const labelList& faceAnchorLevel,
                const labelList& edgeMidPoint,
                const label cellI,
                polyTopoChange& ref
            ) const;


        // Topological change helper functions

            //- Get cell added to point of cellI (if any)
            label getAnchorCell
            (
                const labelListList& cellAnchorPoints,
                const labelListList& cellAddedCells,
                const label cellI,
                const label faceI,
                const label pointI
            ) const;

            //- Set new owner and neighbour (in unspecified order) of pointI
            //  on faceI
            void setNewFaceNeighbours
            (
                const labelListList& cellAnchorPoints,
                const labelListList& cellAddedCells,
                const label faceI,
                const label pointI,

                label& own,
                label& nei
            ) const;

            //- Find index of point with wantedLevel, starting from fp
            label findLevel
            (
                const face& f,
                const label startFp,
                const bool searchForward,
                const label wantedLevel
            ) const;

            //- Store in maps correspondence from midpoint to anchors and
            //  faces. Used when creating internal faces
            label storeMidPointInfo
            (
                const labelListList& cellAnchorPoints,
                const labelListList& cellAddedCells,
                const labelList& cellMidPoint,
                const labelList& edgeMidPoint,
                const label cellI,
                const label faceI,
                const bool faceOrder,
                const label midPointI,
                const label anchorPointI,
                const label faceMidPointI,

                Map<edge>& midPointToAnchors,
                Map<edge>& midPointToFaceMids,
                polyTopoChange& ref
            ) const;

            //- If p0 and p1 are existing vertices check if edge is split and insert
            //  splitPoint. Used with storing mid point
            void insertEdgeSplit
            (
                const labelList& edgeMidPoint,
                const label p0,
                const label p1,
                DynamicList<label>& verts
            ) const;


        // Copy control

            //- Disallow default bitwise copy construct
            polyhedralRefinement(const polyhedralRefinement&);

            //- Disallow default bitwise assignment
            void operator=(const polyhedralRefinement&);


protected:

    // Protected Pure Virtual Member Functions

        // Global topology modification functions (operate on whole polyMesh)

            //- Set refinement instruction
            virtual void setRefinement
            (
                polyTopoChange& meshMod,
                const labelList& cellsToRefine
            ) const;

            //- Set unrefinement instruction
            virtual void setUnrefinement
            (
                polyTopoChange& meshMod,
                const labelList& splitPointsToUnrefine
            ) const;


public:

    //- Runtime type information
    TypeName("polyhedralRefinement");


    // Constructors

        //- Construct from dictionary
        polyhedralRefinement
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const bool read = true
        );


    //- Destructor
    virtual ~polyhedralRefinement();


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
