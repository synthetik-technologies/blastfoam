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
    Foam::treeNode

Description
    Class to implement octree.

    Holds the pointers to sub-octants. These are either other treeNodes or
    treeLeafs. The treeLeafs hold the actual data as a list of indices into
    octreeData.

Note
    To prevent calculation errors all bounding boxes used in octrees are
    calculated only once.

    The pointers to either treeNode/treeLeaf are implemented 'by hand'
    (explicitly marking type) instead of using a proper virtual mechanism
    to save some space in the treeLeaves.

SourceFiles
    treeNode.C

\*---------------------------------------------------------------------------*/

#ifndef treeNode_H
#define treeNode_H

#include "treeBoundBoxList.H"
#include "treeElem.H"
#include "linePointRef.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// class intersection;
template<class Type> class octree;
template<class Type> class treeLeaf;


// Forward declaration of friend functions and operators
template<class Type> class treeNode;

template<class Type> Istream& operator>>(Istream&, treeNode<Type>&);
template<class Type> Ostream& operator<<(Ostream&, const treeNode<Type>&);


/*---------------------------------------------------------------------------*\
                        Class treeNodeName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(treeNode);


/*---------------------------------------------------------------------------*\
                          Class treeNode Declaration
\*---------------------------------------------------------------------------*/

template <class Type>
class treeNode
:
    public treeElem<Type>,
    public treeNodeName
{
    // Private data

        //- Position of the midpoint
        const point mid_;

        //- Type stored in subNodes_
        unsigned char subNodeTypes_;

        //- Pointers to sub treeNode or treeLeaf
        treeElem<Type>* subNodes_[8];

        //- Constant valid for whole subNode/leaf
        label volType_;

    // Static data members

        //- leaf offset for octant index
        static const label leafOffset;


    // Private Member Functions

        //- mark pointer to subnode as being a treeNode*
        void setAsNode(const label octant);

        //- mark pointer to subnode as being a treeLeaf*
        void setAsLeaf(const label octant);

        //- Set pointer to sub node
        void setNodePtr(const label octant, treeElem<Type>* treeNodePtr);

        //- Set pointer to sub leaf
        void setLeafPtr(const label octant, treeElem<Type>* treeLeafPtr);

        //- Set type of octant
        void setVolType(const label octant, const label type);

        //- Get type of octant
        inline label getVolType(const label octant) const;

        //- Find first leaf on line start-end. Updates start.
        const treeLeaf<Type>* findLeafLineOctant
        (
            const int level,
            const Type& shapes,
            const label octant,
            const vector& direction,
            point& start,
            const point& end
        ) const;


        //- Print spaces
        static void space(Ostream&, const label);

        //- Disallow default bitwise copy construct
        treeNode(const treeNode&);

        //- Disallow default bitwise assignment
        void operator=(const treeNode&);


public:

    // Constructors

        //- Construct from components
        treeNode(const treeBoundBox&);

        //- Construct from Istream
        treeNode(Istream&);


    // Destructor

        ~treeNode();


    // Member Functions

        // Access

            //- The midpoint position
            inline const point& midpoint() const;

            //- array of 8 subNodes/leaves
            inline treeElem<Type>* const* subNodes() const;

            //- octant contains pointer to treeNode(1) or treeLeaf(0)
            inline label isNode(const label octant) const;

            //- Get pointer to sub node
            inline treeNode<Type>* getNodePtr(const label octant) const;

            //- Get pointer to sub leaf
            inline treeLeaf<Type>* getLeafPtr(const label octant) const;

        // Edit

            //- Take list of shapes and distribute over the 8 octants
            void distribute
            (
                const label,
                octree<Type>&,
                const Type& shapes,
                const labelList&
            );

            //- Distribute at certain level only
            void redistribute
            (
                const label,
                octree<Type>&,
                const Type& shapes,
                const label
            );

            //- Set type of subnodes
            label setSubNodeType
            (
                const label level,
                octree<Type>& top,
                const Type& shapes
            );

        // Search

            //- Find type of node sample is in. Used for inside/outside
            //  determination
            label getSampleType
            (
                const label level,
                const octree<Type>& top,
                const Type& shapes,
                const point& sample
            ) const;

            //- Find index of shape containing sample.
            label find
            (
                const Type& shapes,
                const point& sample
            ) const;

            //- Find tightest bounding box around sample which is guaranteed
            //  to hold at least one cell.
            //  Current best bb in tightest,
            //  returns true if newTightest has changed, 0 otherwise.
            bool findTightest
            (
                const Type& shapes,
                const point& sample,
                treeBoundBox& tightest
            ) const;

            //- Find nearest shape to sample
            //  Returns true if found nearer shape and updates
            //  tightest, tightestI, tightestDist
            bool findNearest
            (
                const Type& shapes,
                const point& sample,
                treeBoundBox& tightest,
                label& tightestI,
                scalar& tightestDist
            ) const;

            //- Find nearest shape to line
            //  Returns true if found nearer shape and updates nearest and
            //  tightest
            bool findNearest
            (
                const Type& shapes,
                const linePointRef& ln,
                treeBoundBox& tightest,
                label& tightestI,   // index of nearest shape
                point& linePoint,   // nearest point on line
                point& shapePoint   // nearest point on shape
            ) const;

            //- Find shapes not outside box. Return true if anything found.
            bool findBox
            (
                const Type& shapes,
                const treeBoundBox& bb,
                labelHashSet& elements
            ) const;

            //- Find treeLeaves intersecting line segment [start..end]
            //  Updates: start
            const treeLeaf<Type>* findLeafLine
            (
                const int level,
                const Type& shapes,
                point& start,
                const point& end
            ) const;


            //- Collect all treeLeafs in leafArray. leafIndex points to first
            //  empty slot in leafArray and gets updated.
            void findLeaves
            (
                List<treeLeaf<Type>*>& leafArray,
                label& leafIndex
            ) const;

            //- Same but for const.
            void findLeaves
            (
                List<const treeLeaf<Type>*>& leafArray,
                label& leafIndex
            ) const;


        // Write

            //- Print contents of node.
            void printNode
            (
                Ostream& os,
                const label level
            ) const;

            //- Write subleafs in OBJ format.
            void writeOBJ
            (
                Ostream& os,
                const int level,
                label& vertNo
            ) const;


    // IOstream Operators

        friend Istream& operator>> <Type> (Istream&, treeNode<Type>&);
        friend Ostream& operator<< <Type> (Ostream&, const treeNode<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


#include "treeNodeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "treeNode.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
