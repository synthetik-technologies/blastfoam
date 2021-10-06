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
    Foam::octreeLine

Description
    Iterates over intersections of line with octree leaf elements.

    Used as in
    @code
        octree<octreeDataFace> oc( .. );

        octreeLine<octreeDataFace> lineSearch(oc, pStart, pEnd);

        while (lineSearch.getIntersection())
        {
            const point& pt = lineSearch.hitInfo().hitPoint();
            ..
        }
    @endcode

SourceFiles
    octreeLine.C

\*---------------------------------------------------------------------------*/

#ifndef octreeLine_H
#define octreeLine_H

#include "boolList.H"
#include "point.H"
#include "pointHitSort.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class Type> class octree;
template<class Type> class treeLeaf;


/*---------------------------------------------------------------------------*\
                           Class octreeLine Declaration
\*---------------------------------------------------------------------------*/

template <class Type>
class octreeLine
{
    // Private data

        //- Octree reference
        const octree<Type>& tree_;

        //- Start of segment
        const point startPoint_;

        //- End of segment
        const point endPoint_;

        //- Start moved into bb
        point realStartPoint_;

        //- Exit point of intersection with current treeLeaf
        point leafExitPoint_;

        //- Current treeLeaf to be searched in.
        const treeLeaf<Type>* currentLeaf_;

        //- Sorted list of intersections
        List<pointHitSort> sortedIntersections_;

        //- index of last hit in previous treeLeaf. Used so if shape double
        //  it does not get counted twice. Note is not ok for concave shapes
        label lastElem_;

        //- Current hit: index in sortedIntersections_
        label sortedI_;

    // Private Member Functions

        //- Calculate sorted list of intersections
        void calcSortedIntersections();

        //- Searches for leaf with intersected elements.
        //  Return true if found; false otherwise.
        //  Sets currentLeaf_ and sortedIntersections_
        bool getNextLeaf();

public:

    // Constructors

        //- Construct from components
        octreeLine
        (
            const octree<Type>& tree,
            const point& startPoint,
            const point& endPoint
        );


    // Destructor

        ~octreeLine();


    // Member Functions

        const octree<Type>& tree() const
        {
            return tree_;
        }

        const point& leafExitPoint() const
        {
            return leafExitPoint_;
        }

        const point& endPoint() const
        {
            return endPoint_;
        }

        const point& startPoint() const
        {
            return startPoint_;
        }

        const treeLeaf<Type>* currentLeaf() const
        {
            return currentLeaf_;
        }

        const List<pointHitSort>& sortedIntersections() const
        {
            return sortedIntersections_;
        }

        label hitIndex() const
        {
            return sortedIntersections_[sortedI_].index();
        }

        const pointHit& hitInfo() const
        {
            return sortedIntersections_[sortedI_].inter();
        }


        //- go to next intersection. Return false if no intersections.
        bool getIntersection();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "octreeLine.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
