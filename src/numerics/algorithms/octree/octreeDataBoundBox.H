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
    Foam::octreeDataBoundBox

Description
    Holds data for octree to work on simple bounding boxes

    For example, calculate (in calcNearest) the correct intersection point
    with a bounding box.

    Using bounding boxes simplify the usage of octree when only
    primitivePatches entities are available (eg: from inside ggi interpolators)

SourceFiles
    octreeDataBoundBox.C

Author
    Martin Beaudoin, Hydro-Quebec, (2010)

\*---------------------------------------------------------------------------*/

#ifndef octreeDataBoundBox_H
#define octreeDataBoundBox_H

#include "treeBoundBoxList.H"
#include "faceList.H"
#include "point.H"
#include "className.H"
#include "linePointRef.H"
// #include "tolerancesSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class primitiveMesh;
template<class Type> class octree;
class polyPatch;

/*---------------------------------------------------------------------------*\
                           Class octreeDataBoundBox Declaration
\*---------------------------------------------------------------------------*/

class octreeDataBoundBox
{
    // Static data

        //- Tolerance on linear dimensions
        static const /*debug::tolerancesSwitch*/ scalar tol;


    // Private data

        //- Bounding boxes
        treeBoundBoxList allBb_;


    // Private Member Functions


public:

    // Declare name of the class and its debug switch
    ClassName("octreeDataBoundBox");


    // Constructors

        //- Construct from list of precomputed of bound boxes
        explicit octreeDataBoundBox
        (
            const treeBoundBoxList& bbL
        );

        //- Construct as copy
        octreeDataBoundBox(const octreeDataBoundBox&);


    // Destructor

        ~octreeDataBoundBox();


    // Member Functions

        // Access

            const treeBoundBoxList& allBb() const
            {
                return allBb_;
            }

            label size() const
            {
                return allBb_.size();
            }


        // Search

            //- Get type of sample
            label getSampleType
            (
                const octree<octreeDataBoundBox>&,
                const point&
            ) const;

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
            bool intersects
            (
                const label index,
                const point& start,
                const point& end,
                point& intersectionPoint
            ) const;

            //- Sets newTightest to bounding box (and returns true) if
            //  nearer to sample than tightest bounding box. Otherwise
            //  returns false.
            bool findTightest
            (
                const label index,
                const point& sample,
                treeBoundBox& tightest
            ) const;

            //- Given index get unit normal and calculate (numerical) sign
            //  of sample.
            //  Used to determine accuracy of calcNearest or inside/outside.
            scalar calcSign
            (
                const label index,
                const point& sample,
                vector& n
            ) const;

            //- Calculates nearest (to sample) point in shape.
            //  Returns point and mag(nearest - sample). Returns GREAT if
            //  sample does not project onto (triangle decomposition) of face.
            scalar calcNearest
            (
                const label index,
                const point& sample,
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

            //- Write shape at index
            void write(Ostream& os, const label index) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
