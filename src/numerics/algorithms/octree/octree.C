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

Description

\*---------------------------------------------------------------------------*/

#include "octree.H"
#include "treeLeaf.H"
#include "treeNode.H"
#include "cpuTime.H"
#include "linePointRef.H"
#include "pointIndexHit.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template <class Type>
Foam::string Foam::octree<Type>::volType(const label type)
{
    if (type == UNKNOWN)
    {
        return "unknown";
    }
    else if (type == MIXED)
    {
        return "mixed";
    }
    else if (type == INSIDE)
    {
        return "inside";
    }
    else if (type == OUTSIDE)
    {
        return "outside";
    }
    else
    {
        FatalErrorIn("volType(const label)") << "type:" << type
            << " unknown." << abort(FatalError);

        return "dummy";
    }
}


// Determine inside/outside status of vector compared to geometry-based normal
template <class Type>
Foam::label Foam::octree<Type>::getVolType
(
    const vector& geomNormal,
    const vector& vec
)
{
    scalar sign = geomNormal & vec;

    if (sign >= 0)
    {
        return octree<Type>::OUTSIDE;
    }
    else
    {
        return octree<Type>::INSIDE;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class Type>
Foam::octree<Type>::octree
(
    const treeBoundBox& octreeBb,
    const Type& shapes,
    const label minNLevels,
    const scalar maxLeafRatio,
    const scalar maxShapeRatio
)
:
    topNode_(new treeNode<Type>(octreeBb)),
    shapes_(shapes),
    octreeBb_(octreeBb),
    maxLeafRatio_(maxLeafRatio),
    maxShapeRatio_(maxShapeRatio),
    minNLevels_(minNLevels),
    deepestLevel_(0),
    nEntries_(0),
    nNodes_(0),
    nLeaves_(0),
    endIter_(*this, -1),
    endConstIter_(*this, -1)
{
    cpuTime timer;

    setNodes(nNodes() + 1);

    const label nShapes = shapes_.size();

    labelList indices(nShapes);
    for (label i = 0; i < nShapes; i++)
    {
        indices[i] = i;
    }

    // Create initial level (0) of subLeaves
    if (debug & 1)
    {
        Pout<< "octree : --- Start of Level " << deepestLevel_
            << " ----" << endl;
    }
    topNode_->distribute(0, *this, shapes_, indices);

    if (debug & 1)
    {
        printStats(Pout);
        Pout<< "octree : --- End of Level " << deepestLevel_
            << " ----" << endl;
    }

    // Breadth first creation of tree
    // Stop if: - level above minlevel and
    //                - less than so many cells per endpoint
    //                  (so bottom level is fine enough)
    //                - every shape mentioned in only so many
    //                  leaves. (if shape bb quite big it will end up
    //                  in lots of leaves).
    //          - no change in number of leaves
    //            (happens if leafs fine enough)
    // This guarantees that tree
    //  - is fine enough (if minLevel > 0)
    //  - has some guaranteed maximum size (maxShapeRatio)

    label oldNLeaves = -1;  // make test below pass first time.
    deepestLevel_ = 1;
    while
    (
        (deepestLevel_ <= minNLevels_)
     || (
            (nEntries() > maxLeafRatio * nLeaves())    // shapes per leaf
         && (nEntries() < maxShapeRatio * nShapes)     // entries per shape
        )
    )
    {
        if (deepestLevel_ >= maxNLevels)
        {
            if (debug & 1)
            {
                Pout<< "octree : exiting since maxNLevels "
                    << maxNLevels << " reached" << endl;
            }
            break;
        }

        if (oldNLeaves == nLeaves())
        {
            if (debug & 1)
            {
                Pout<< "octree : exiting since nLeaves does not change"
                    << endl;
            }
            break;
        }
        if (debug & 1)
        {
            Pout<< "octree : --- Start of Level " << deepestLevel_
                << " ----" << endl;
        }

        oldNLeaves = nLeaves();

        topNode_->redistribute
        (
            1,
            *this,
            shapes_,
            deepestLevel_
        );

        if (debug & 1)
        {
            printStats(Pout);

            Pout<< "octree : --- End of Level " << deepestLevel_
                << " ----" << endl;
        }

        deepestLevel_++;
    }


    if (debug & 1)
    {
        Pout<< "octree : Constructed octree in = "
        << timer.cpuTimeIncrement()
        << " s\n" << endl << endl;
    }

    // Set volume type of non-treeleaf nodes.
    topNode_->setSubNodeType(0, *this, shapes_);

    if (debug & 1)
    {
        Pout<< "octree : Added node information to octree in  = "
        << timer.cpuTimeIncrement()
        << " s\n" << endl << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class Type>
Foam::octree<Type>::~octree()
{
    delete topNode_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class Type>
Foam::label Foam::octree<Type>::getSampleType(const point& sample) const
{
    return topNode_->getSampleType(0, *this, shapes_, sample);
}


template <class Type>
Foam::label Foam::octree<Type>::find(const point& sample) const
{
    return topNode_->find(shapes_, sample);
}


template <class Type>
bool Foam::octree<Type>::findTightest
(
    const point& sample,
    treeBoundBox& tightest
) const
{
    return topNode_->findTightest
    (
        shapes_,
        sample,
        tightest
    );
}


template <class Type>
Foam::label Foam::octree<Type>::findNearest
(
    const point& sample,
    treeBoundBox& tightest,
    scalar& tightestDist
) const
{
    label tightestI = -1;

    if (debug & 4)
    {
        Pout<< "octree::findNearest : searching for nearest for "
            << "sample:" << sample
            << "  tightest:" << tightest << endl;
    }

    topNode_->findNearest
    (
        shapes_,
        sample,
        tightest,
        tightestI,
        tightestDist
    );

    if (debug & 4)
    {
        Pout<< "octree::findNearest : found nearest for "
            << "sample:" << sample << " with "
            << "  tightestI:" << tightestI
            << "  tightest:" << tightest
            << "  tightestDist:" << tightestDist
            << endl;
    }

    return tightestI;
}


template <class Type>
Foam::label Foam::octree<Type>::findNearest
(
    const linePointRef& ln,
    treeBoundBox& tightest,
    point& linePoint,
    point& shapePoint
) const
{
    // Start off from miss with points at large distance apart.
    label tightestI = -1;
    linePoint = point(-GREAT, -GREAT, -GREAT);
    shapePoint = point(GREAT, GREAT, GREAT);

    topNode_->findNearest
    (
        shapes_,
        ln,
        tightest,
        tightestI,
        linePoint,
        shapePoint
    );

    return tightestI;
}


template <class Type>
Foam::labelList Foam::octree<Type>::findBox(const treeBoundBox& bb) const
{
    // Storage for labels of shapes inside bb. Size estimate.
    labelHashSet elements(100);

    topNode_->findBox(shapes_, bb, elements);

    return elements.toc();
}


template <class Type>
Foam::pointIndexHit Foam::octree<Type>::findLine
(
    const point& treeStart,
    const point& treeEnd
) const
{
    // Initialize to a miss
    pointIndexHit hitInfo(false, treeStart, -1);

    const vector dir(treeEnd - treeStart);

    // Current line segment to search
    point start(treeStart);
    point end(treeEnd);

    while (true)
    {
        // Find nearest treeLeaf intersected by line
        point leafIntPoint;

        const treeLeaf<Type>* leafPtr = findLeafLine
        (
            start,
            end,
            leafIntPoint
        );

        if (!leafPtr)
        {
            // Reached end of string of treeLeaves to be searched. Return
            // whatever we have found so far.
            break;
        }

        // Inside treeLeaf find nearest intersection
        scalar minS = GREAT;

        const labelList& indices = leafPtr->indices();

        forAll(indices, elemI)
        {
            label index = indices[elemI];

            point pt;
            bool hit = shapes().intersects(index, start, end, pt);

            if (hit)
            {
                // Check whether intersection nearer p
                scalar s = (pt - treeStart) & dir;

                if (s < minS)
                {
                    minS = s;

                    // Update hit info
                    hitInfo.setHit();
                    hitInfo.setPoint(pt);
                    hitInfo.setIndex(index);

                    // Update segment to search
                    end = pt;
                }
            }
        }

        if (hitInfo.hit())
        {
            // Found intersected shape.
            break;
        }

        // Start from end of current leaf
        start = leafIntPoint;
    }

    return hitInfo;
}


template <class Type>
Foam::pointIndexHit Foam::octree<Type>::findLineAny
(
    const point& start,
    const point& end
) const
{
    // Initialize to a miss
    pointIndexHit hitInfo(false, start, -1);

    // Start of segment in current treeNode.
    point p(start);
    while (true)
    {
        // Find treeLeaf intersected by line
        point leafIntPoint;

        const treeLeaf<Type>* leafPtr = findLeafLine(p, end, leafIntPoint);

        if (!leafPtr)
        {
            // Reached end of string of treeLeaves to be searched. Return
            // whatever we have found so far.
            break;
        }

        // Inside treeLeaf find any intersection

        const labelList& indices = leafPtr->indices();

        forAll(indices, elemI)
        {
            label index = indices[elemI];

            point pt;
            bool hit = shapes().intersects
            (
                index,
                p,
                end,
                pt
            );

            if (hit)
            {
                hitInfo.setHit();
                hitInfo.setPoint(pt);
                hitInfo.setIndex(index);

                break;
            }
        }

        if (hitInfo.hit())
        {
            // Found intersected shape.
            break;
        }

        // Start from end of current leaf
        p = leafIntPoint;
    }

    return hitInfo;
}


template <class Type>
const Foam::treeLeaf<Type>* Foam::octree<Type>::findLeafLine
(
    const point& start,
    const point& end,
    point& leafIntPoint
) const
{
    // returns first found cube along line

    if (debug & 2)
    {
        Pout<< "octree::findLeafLine : searching for shapes on line "
            << "start:" << start
            << "  end:" << end << endl;
    }

    // If start is outside project onto top cube
    if (octreeBb_.contains(start))
    {
        leafIntPoint = start;
    }
    else
    {
        if (!octreeBb_.intersects(start, end, leafIntPoint))
        {
            if (debug & 2)
            {
                Pout<< "octree::findLeafLine : start outside domain but does"
                    << " not intersect : start:"
                    << start << endl;
            }
            return nullptr;
        }

        if (debug & 2)
        {
            Pout<< "octree::findLeafLine : start propagated to inside"
                   " domain : "
                << leafIntPoint << endl;
        }
    }

    // Normal action: find next intersection along line
    const treeLeaf<Type>* leafPtr = topNode_->findLeafLine
    (
        0,
        shapes_,
        leafIntPoint,
        end
    );

    if (debug & 2)
    {
        Pout<< "returning from octree::findLeafLine with "
            << "leafIntersection:" << leafIntPoint
            << "  leafPtr:" << long(leafPtr) << endl;
    }

    return leafPtr;
}


template <class Type>
void Foam::octree<Type>::writeOBJ
(
    Ostream& os,
    label& vertNo
) const
{
    scalar minx = octreeBb_.min().x();
    scalar miny = octreeBb_.min().y();
    scalar minz = octreeBb_.min().z();

    scalar maxx = octreeBb_.max().x();
    scalar maxy = octreeBb_.max().y();
    scalar maxz = octreeBb_.max().z();

    os << "v " << minx << " " << miny << " " << minz << endl;
    os << "v " << maxx << " " << miny << " " << minz << endl;
    os << "v " << maxx << " " << maxy << " " << minz << endl;
    os << "v " << minx << " " << maxy << " " << minz << endl;

    os << "v " << minx << " " << miny << " " << maxz << endl;
    os << "v " << maxx << " " << miny << " " << maxz << endl;
    os << "v " << maxx << " " << maxy << " " << maxz << endl;
    os << "v " << minx << " " << maxy << " " << maxz << endl;

    // Bottom face
    os << "l " << vertNo + 1 << " " << vertNo + 2 << endl;
    os << "l " << vertNo + 2 << " " << vertNo + 3 << endl;
    os << "l " << vertNo + 3 << " " << vertNo + 4 << endl;
    os << "l " << vertNo + 4 << " " << vertNo + 1 << endl;

    // Top face
    os << "l " << vertNo + 5 << " " << vertNo + 6 << endl;
    os << "l " << vertNo + 6 << " " << vertNo + 7 << endl;
    os << "l " << vertNo + 7 << " " << vertNo + 8 << endl;
    os << "l " << vertNo + 8 << " " << vertNo + 5 << endl;

    // Edges from bottom to top face
    os << "l " << vertNo + 1 << " " << vertNo + 5 << endl;
    os << "l " << vertNo + 2 << " " << vertNo + 6 << endl;
    os << "l " << vertNo + 3 << " " << vertNo + 7 << endl;
    os << "l " << vertNo + 4 << " " << vertNo + 8 << endl;

    vertNo += 8;

    topNode_->writeOBJ(os, 1, vertNo);
}


template <class Type>
void Foam::octree<Type>::printStats(Ostream& os) const
{
    os  << "Statistics after iteration " << deepestLevel() << ':' << endl
        << "  nShapes  :" << shapes().size() << endl
        << "  nNodes   :" << nNodes() << endl
        << "  nLeaves  :" << nLeaves() << endl
        << "  nEntries :" << nEntries() << endl;

    if (nLeaves() && shapes().size())
    {
        os
            << "  Cells per leaf :"
            << scalar(nEntries())/nLeaves()
            << nl
            << "  Every cell in  :"
            << scalar(nEntries())/shapes().size() << " cubes"
            << endl;
    }
}


// * * * * * * * * * * * * * * * STL iterator  * * * * * * * * * * * * * * * //

// Construct from a octree. Set index at end.
template <class Type>
Foam::octree<Type>::iterator::iterator(octree<Type>& oc)
:
    octree_(oc),
    curLeaf_(oc.nLeaves())
{
    leaves_.setSize(0);
}


// Construct from octree. Set index.
template <class Type>
Foam::octree<Type>::iterator::iterator(octree<Type>& oc, label index)
:
    octree_(oc),
    curLeaf_(index)
{
    if (index >= 0)
    {
        leaves_.setSize(oc.nLeaves());

        label leafIndex = 0;
        octree_.topNode()->findLeaves(leaves_, leafIndex);

        if (leafIndex != oc.nLeaves())
        {
            FatalErrorIn
            (
                "octree::iterator::iterator"
                "(octree&, label)"
            )
            << "Traversal of tree returns : " << leafIndex << endl
            << "Statistics of octree say  : " << oc.nLeaves() << endl
            << abort(FatalError);
        }
    }
}


template <class Type>
void Foam::octree<Type>::iterator::operator=(const iterator& iter)
{
    if ((curLeaf_ < 0) && (iter.curLeaf_ >= 0))
    {
        FatalErrorIn
        (
            "octree::iterator::operator="
            "(const iterator&)"
        )
        << "lhs : " << curLeaf_ << endl
        << "rhs : " << iter.curLeaf_ << endl
        << abort(FatalError);
    }
    curLeaf_ = iter.curLeaf_;
}


template <class Type>
bool Foam::octree<Type>::iterator::operator==(const iterator& iter) const
{
    label index1 =
        (curLeaf_ >= 0 ? curLeaf_ : octree_.nLeaves());
    label index2 =
        (iter.curLeaf_ >= 0 ? iter.curLeaf_ : iter.octree_.nLeaves());

    return index1 == index2;
}


template <class Type>
bool Foam::octree<Type>::iterator::operator!=(const iterator& iter) const
{
    return !(iterator::operator==(iter));
}


template <class Type>
Foam::treeLeaf<Type>& Foam::octree<Type>::iterator::operator*()
{
    return *leaves_[curLeaf_];
}


template <class Type>
typename Foam::octree<Type>::iterator&
Foam::octree<Type>::iterator::operator++()
{
    curLeaf_++;
    return *this;
}


template <class Type>
typename Foam::octree<Type>::iterator
Foam::octree<Type>::iterator::operator++(int)
{
    iterator tmp = *this;
    ++*this;
    return tmp;
}


template <class Type>
typename Foam::octree<Type>::iterator
Foam::octree<Type>::begin()
{
    return iterator(*this, 0);
}


template <class Type>
const typename Foam::octree<Type>::iterator&
Foam::octree<Type>::end()
{
    return octree<Type>::endIter_;
}


// * * * * * * * * * * * * * * STL const_iterator  * * * * * * * * * * * * * //

// Construct for a given octree
template <class Type>
Foam::octree<Type>::const_iterator::const_iterator(const octree<Type>& oc)
:
    octree_(oc),
    curLeaf_(oc.nLeaves())
{
    leaves_.setSize(oc.nLeaves());
}


// Construct for a given octree
template <class Type>
Foam::octree<Type>::const_iterator::const_iterator
(
    const octree<Type>& oc,
    label index
)
:
    octree_(oc),
    curLeaf_(index)
{
    if (index >= 0)
    {
        leaves_.setSize(oc.nLeaves());

        label leafIndex = 0;
        octree_.topNode()->findLeaves(leaves_, leafIndex);

        if (leafIndex != oc.nLeaves())
        {
            FatalErrorIn
            (
                "octree::const_iterator::const_iterator"
                "(octree&, label)"
            )
            << "Traversal of tree returns : " << leafIndex << endl
            << "Statistics of octree say  : " << oc.nLeaves() << endl
            << abort(FatalError);
        }
    }
}


template <class Type>
void Foam::octree<Type>::const_iterator::operator=(const const_iterator& iter)
{
    if ((curLeaf_ < 0) && (iter.curLeaf_ >= 0))
    {
        FatalErrorIn
        (
            "octree::const_iterator::operator="
            "(const const_iterator&)"
        )
        << "lhs : " << curLeaf_ << endl
        << "rhs : " << iter.curLeaf_ << endl
        << abort(FatalError);
    }
    curLeaf_ = iter.curLeaf_;
    curLeaf_ = iter.curLeaf_;
}


template <class Type>
bool Foam::octree<Type>::const_iterator::operator==
(
    const const_iterator& iter
) const
{
    label index1 =
        (curLeaf_ >= 0 ? curLeaf_ : octree_.nLeaves());
    label index2 =
        (iter.curLeaf_ >= 0 ? iter.curLeaf_ : iter.octree_.nLeaves());

    return index1 == index2;
}


template <class Type>
bool Foam::octree<Type>::const_iterator::operator!=
(
    const const_iterator& iter
) const
{
    return !(const_iterator::operator==(iter));
}


template <class Type>
const Foam::treeLeaf<Type>& Foam::octree<Type>::const_iterator::operator*()
{
    return *leaves_[curLeaf_];
}


template <class Type>
typename Foam::octree<Type>::const_iterator&
Foam::octree<Type>::const_iterator::operator++()
{
    curLeaf_++;
    return *this;
}


template <class Type>
typename Foam::octree<Type>::const_iterator
Foam::octree<Type>::const_iterator::operator++(int)
{
    const_iterator tmp = *this;
    ++*this;
    return tmp;
}


template <class Type>
typename Foam::octree<Type>::const_iterator
Foam::octree<Type>::begin() const
{
    return const_iterator(*this, 0);
}


template <class Type>
typename Foam::octree<Type>::const_iterator
Foam::octree<Type>::cbegin() const
{
    return const_iterator(*this, 0);
}


template <class Type>
const typename Foam::octree<Type>::const_iterator&
Foam::octree<Type>::end() const
{
    return octree<Type>::endConstIter_;
}


template <class Type>
const typename Foam::octree<Type>::const_iterator&
Foam::octree<Type>::cend() const
{
    return octree<Type>::endConstIter_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template <class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const octree<Type>& oc)
{
    return os << token::BEGIN_LIST
        //<< token::SPACE << oc.shapes_
        << token::SPACE << oc.octreeBb_
        << token::SPACE << oc.maxLeafRatio_
        << token::SPACE << oc.maxShapeRatio_
        << token::SPACE << oc.minNLevels_
        << token::SPACE << oc.deepestLevel_
        << token::SPACE << oc.nEntries_
        << token::SPACE << oc.nNodes_
        << token::SPACE << oc.nLeaves_
        << token::SPACE << *oc.topNode_
        << token::SPACE << token::END_LIST;
}


// ************************************************************************* //
