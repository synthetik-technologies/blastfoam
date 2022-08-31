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
    Foam::octreeDataCell

Description
    Encapsulation of data needed to search in/for cells.

    Used to find the cell containing a point (e.g. cell-cell mapping).

SourceFiles
    octreeDataCell.C

\*---------------------------------------------------------------------------*/

#ifndef octreeDataCell_H
#define octreeDataCell_H

#include "treeBoundBoxList.H"
#include "labelList.H"
#include "linePointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
template<class Type> class octree;

/*---------------------------------------------------------------------------*\
                       Class octreeDataCell Declaration
\*---------------------------------------------------------------------------*/

class octreeDataCell
{
    // Private data

        const polyMesh& mesh_;

        labelList cellLabels_;

        treeBoundBoxList bbs_;


public:

    // Constructors

        //- Construct from components.
        octreeDataCell
        (
            const polyMesh&,
            const labelList& cellLabels,
            const treeBoundBoxList& bbs
        );

        //- Construct from mesh. Uses all cells in mesh.
        octreeDataCell(const polyMesh&);


    // Member Functions

        // Access

            const labelList& cellLabels() const
            {
                return cellLabels_;
            }

            const polyMesh& mesh() const
            {
                return mesh_;
            }

            const treeBoundBoxList& allBb() const
            {
                return bbs_;
            }

            label size() const
            {
                return bbs_.size();
            }

        // Search

            //- Get type of sample
            label getSampleType(octree<octreeDataCell>&, const point&) const;

            //- Does (bb of) shape at index overlap bb
            bool overlaps
            (
                const label index,
                const treeBoundBox& sampleBb
            ) const;

            //- Does shape at index contain sample
            bool contains
            (
                const label index,
                const point& sample
            ) const;

            //- Segment (from start to end) intersection with shape
            //  at index. If intersects returns true and sets intersectionPoint
            // BUG: not yet done.
            bool intersects
            (
                const label index,
                const point& start,
                const point& end,
                point& intersectionPoint
            ) const;

            //- Sets newTightest to bounding box (and returns true) if
            //  nearer to sample than tightest bounding box. Otherwise
            //  returns false
            bool findTightest
            (
                const label index,
                const point& sample,
                treeBoundBox& tightest
            ) const;

            //- Given index get unit normal and calculate (numerical) sign
            //  of sample.
            //  Used to determine accuracy of calcNearest or inside/outside.
            //  Note: always returns GREAT since no inside/outside.
            scalar calcSign
            (
                const label index,
                const point& sample,
                vector& n
            ) const;

            //- Calculates nearest (to sample) point in shape.
            //  Returns point and mag(nearest - sample)
            scalar calcNearest
            (
                const Foam::label index,
                const Foam::point& sample,
                point& nearest
            ) const;

            //- Calculates nearest (to line segment) point in shape.
            //  Returns distance and both point.
            scalar calcNearest
            (
                const label index,
                const linePointRef& ln,
                point& linePt,          // nearest point on line
                point& shapePt          // nearest point on shape
            ) const;



        // Write

            // Write shape at index
            void write(Ostream& os, const label index) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
