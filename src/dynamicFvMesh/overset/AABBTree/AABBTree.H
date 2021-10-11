/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::AABBTree

Description
    Templated tree of axis-aligned bounding boxes (AABB)

    Designed to be templated on either faces or cells, the AABBTree will
    decompose the input into a tree of AABB's.  The maximum number of tree
    levels and minimum number of objects per leaf are provided on construction,
    and the contents (addressing) is stored.

SourceFiles
    AABBTree.C

\*---------------------------------------------------------------------------*/

#ifndef AABBTree_H
#define AABBTree_H

#include "labelList.H"
#include "labelPair.H"
#include "DynamicList.H"
#include "pointField.H"
#include "treeBoundBox.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type>
class AABBTree;

template<class Type>
Istream& operator>>(Istream&, AABBTree<Type>&);

template<class Type>
Ostream& operator<<(Ostream&, const AABBTree<Type>&);

/*---------------------------------------------------------------------------*\
                          Class AABBTree Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class AABBTree
{

protected:

    // Protected Data

        //- Tolerance
        static scalar tolerance_;

        //- Maximum tree level
        label maxLevel_;

        //- Minimum points per leaf
        label minLeafSize_;

        //- Bounding boxes making up the tree
        List<treeBoundBox> boundBoxes_;

        //- Leaf addressing
        List<labelList> addressing_;


    // Protected Member Functions

        //- Write OBJ file of bounding box
        void writeOBJ
        (
            const bool writeLinesOnly,
            const treeBoundBox& bb,
            label& vertI,
            Ostream& os
        ) const;

        //- Write OBJ for all bounding boxes
        void writeOBJ
        (
            const bool leavesOnly,
            const bool writeLinesOnly,
            const treeBoundBox& bb,
            const label nodeI,
            const List<Pair<treeBoundBox>>& bbs,
            const List<Pair<label>>& nodes,
            label& vertI,
            Ostream& os
        ) const;

        //- Create the bounding boxes by interrogating points
        void createBoxes
        (
            const bool equalBinSize,
            const label level,
            const List<Type>& objects,
            const pointField& points,
            const DynamicList<label>& objectIDs,
            const treeBoundBox& bb,
            const label nodeI,

            DynamicList<Pair<treeBoundBox>>& bbs,
            DynamicList<labelPair>& nodes,
            DynamicList<labelList>& addressing
        ) const;


public:

    // Constructors

        //- Null constructor
        AABBTree();

        //- Construct from components
        //  equalBinSize: divide into equal number of elements or equal span
        AABBTree
        (
            const UList<Type>& objects,
            const pointField& points,
            const bool equalBinSize = true,
            const label maxLevel = 3,
            const label minBinSize = 100
        );


    // Public Member Functions

        // Access

            //- Return the bounding boxes making up the tree
            const List<treeBoundBox>& boundBoxes() const;

            //- Return the contents addressing
            const List<labelList>& addressing() const;


        // Evaluation

            //- Determine whether a point is inside the bounding boxes
            bool pointInside(const point& pt) const;

            //- Determine whether a bounding box overlaps the tree bounding
            //- boxes
            bool overlaps(const boundBox& bbIn) const;


    // IOstream operators

        friend Istream& operator>> <Type>(Istream&, AABBTree&);
        friend Ostream& operator<< <Type>(Ostream&, const AABBTree&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "AABBTree.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
