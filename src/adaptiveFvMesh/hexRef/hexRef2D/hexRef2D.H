/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
21-10-2019  Jeff Heylmun:   Made runtimeSelectable
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::hexRef2D

SourceFiles
    hexRef2D.C

Author
    Oleg Sutyrin
    using hexRef2D code for foam-extend 3.2 by Luca Cornolti

Description
    Refinement of (split) hexes into 4 using directTopoChange based
    on Foam::hexRef8 class

Comment (Luca Cornolti)
    Functions called for un-refining follow the same logic of the related 3D
    case, but they employ the list of central edges instead of cell mid points,
    which do not exist in 2D. For this reason, this part of the code is well
    implemented and could be included in a generalized Foam::hexRef8 class with
    little effort.

    Functions called for refining were implemented following Ahmad
    Baniabedalruhman PhD Thesis with some significant modification.

\*---------------------------------------------------------------------------*/

#ifndef hexRef2D_H
#define hexRef2D_H

#include "hexRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class hexRef2D Declaration
\*---------------------------------------------------------------------------*/

class hexRef2D
:
    public hexRef
{
    // Private Member Functions

        //- Get cell added to point of celli (if any)
        virtual label getAnchorCell
        (
            const labelListList& cellAnchorPoints,
            const labelListList& cellAddedCells,
            const label celli,
            const label facei,
            const label pointi
        ) const;

        //- Store in maps correspondence from midpoint to anchors and faces.
        virtual label storeMidPointInfo
        (
            const labelListList& cellAnchorPoints,
            const labelListList& cellAddedCells,
            const labelList& cellMidPoint,
            const labelList& faceMidPoint,
            const labelList& edgeMidPoint,
            const label celli,
            const label facei,
            const bool faceOrder,
            const label midPointi,
            const label anchorPointi,
            const label faceMidPointi,

            Map<edge>& midPointToAnchors,
            Map<edge>& midPointToFaceMids,
            polyTopoChange& meshMod
        ) const;

        //- Create all internal faces to split celli into 4.
        virtual void createInternalFaces
        (
            const labelListList& cellAnchorPoints,
            const labelListList& cellAddedCells,
            const labelList& cellMidPoint,
            const labelList& faceMidPoint,
            const labelList& faceAnchorLevel,
            const labelList& edgeMidPoint,
            const label celli,
            polyTopoChange& meshMod
        ) const;


public:

    //- Runtime type information
    TypeName("hexRef2D");


    // Constructors

        //- Construct from mesh, read_if_present refinement data
        //  (from write below). If readHistory is true does read_if_present
        //  of refinement history. If false clears all history
        hexRef2D(const polyMesh& mesh, const bool readHistory = true);

        //- Construct from mesh and un/refinement data and optional size of
        //  starting cells
        hexRef2D
        (
            const polyMesh& mesh,
            const labelList& cellLevel,
            const labelList& pointLevel,
            const refinementHistory& history,
            const scalar level0Edge = -great
        );

        //- Construct from mesh and refinement data and optional size of
        //  starting cells
        hexRef2D
        (
            const polyMesh& mesh,
            const labelList& cellLevel,
            const labelList& pointLevel,
            const scalar level0Edge = -great
        );

    // Destructor
        virtual ~hexRef2D();


    // Member Functions

        virtual label extraCells() const
        {
            return 3;
        }

        virtual bool useEdges() const
        {
            return true;
        }

        // Refinement

            //- Insert refinement. All selected cells will be split into 4.
            //  Returns per element in cells the 4 cells they were split into.
            //  Guarantees that the 0th element is the original cell label.
            //  Mapping:
            //  -split cells: 7 new ones get added from original
            //  -split faces: original gets modified; 3 new ones get added
            //               from original
            //  -added internal faces: added from original cell face(if
            //   that was internal) or created out-of-nothing (so will not
            //   get mapped!). Note: could make this inflate from point but
            //   that will allocate interpolation.
            //  -points added to split edge: added from edge start()
            //  -midpoints added: added from cellPoints[0].
            virtual labelListList setRefinement
            (
                const labelList& cells,
                polyTopoChange&
            );


        // Unrefinement (undoing refinement, not arbitrary coarsening)

            //- Return the edges at the centre of top-level split cells
            //  that can be unsplit.
            virtual labelList getSplitPointsEdges() const;

            //- Given proposed
            //  splitEdges to unrefine according to calculate any clashes
            //  (due to 2:1) and return ok list of edges to unrefine.
            //  Either adds points to refine to set (maxSet = true) or
            //  removes points to refine (maxSet = false)
            virtual labelList consistentUnrefinement
            (
                const labelList& edgesToUnrefine,
                const bool maxSet
            ) const;

            //- Remove some refinement. Needs to be supplied output of
            //  consistentUnrefinement. Only call if undoable set.
            //  All 4 edgeCells of a split edge will be combined into
            //  the lowest numbered cell of those 4.
            virtual void setUnrefinement
            (
                const labelList& splitEdgeLabels,
                polyTopoChange&
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
