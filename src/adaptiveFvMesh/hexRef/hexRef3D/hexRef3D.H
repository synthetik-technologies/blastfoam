/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::hexRef3D

Description
    Refinement of (split) hexes using polyTopoChange.

SourceFiles
    hexRef3D.C

\*---------------------------------------------------------------------------*/

#ifndef hexRef3D_H
#define hexRef3D_H

#include "hexRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class hexRef3D Declaration
\*---------------------------------------------------------------------------*/

class hexRef3D
:
    public hexRef
{
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

        //- Create all internal faces to split celli into 8.
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
    TypeName("hexRef3D");


    // Constructors

        //- Construct from mesh, read_if_present refinement data
        //  (from write below). If readHistory is true does read_if_present
        //  of refinement history. If false clears all history
        hexRef3D(const polyMesh& mesh, const bool readHistory = true);

        //- Construct from mesh and un/refinement data and optional size of
        //  starting cells
        hexRef3D
        (
            const polyMesh& mesh,
            const labelList& cellLevel,
            const labelList& pointLevel,
            const refinementHistory& history,
            const scalar level0Edge = -great
        );

        //- Construct from mesh and refinement data and optional size of
        //  starting cells
        hexRef3D
        (
            const polyMesh& mesh,
            const labelList& cellLevel,
            const labelList& pointLevel,
            const scalar level0Edge = -great
        );


    // Destructor
        virtual ~hexRef3D();


    // Member Functions

        // Refinement

            virtual label extraCells() const
            {
                return 7;
            }

            virtual bool useEdges() const
            {
                return false;
            }

            virtual label nDims() const
            {
                return 3;
            }

            //- Insert refinement. All selected cells will be split into 8.
            //  Returns per element in cells the 8 cells they were split into.
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

            //- Return the points at the centre of top-level split cells
            //  that can be unsplit.
            virtual labelList getSplitPointsEdges() const;

            //- Given proposed
            //  splitPoints to unrefine according to calculate any clashes
            //  (due to 2:1) and return ok list of points to unrefine.
            //  Either adds points to refine to set (maxSet = true) or
            //  removes points to refine (maxSet = false)
            virtual labelList consistentUnrefinement
            (
                const labelList& pointsToUnrefine,
                const bool maxSet
            ) const;

            //- Remove some refinement. Needs to be supplied output of
            //  consistentUnrefinement. Only call if undoable set.
            //  All 8 pointCells of a split point will be combined into
            //  the lowest numbered cell of those 8.
            virtual void setUnrefinement
            (
                const labelList& splitPointLabels,
                polyTopoChange&
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
