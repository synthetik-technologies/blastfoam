/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "AABBTree.H"
#include "meshTools.H"
#include "PackedBoolList.H"
//#include "OFstream.H"

template<class Type>
Foam::scalar Foam::AABBTree<Type>::tolerance_ = 1e-4;

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Type>
void Foam::AABBTree<Type>::writeOBJ
(
    const bool writeLinesOnly,
    const treeBoundBox& bb,
    label& vertI,
    Ostream& os
) const
{
    const pointField pts(bb.points());
    for (const point& pt : pts)
    {
        meshTools::writeOBJ(os, pt);
    }

    if (writeLinesOnly)
    {
        for (const edge& e : bb.edges)
        {
            os  << "l " << e[0] + vertI + 1 << ' ' << e[1] + vertI + 1 << nl;
        }
    }
    else
    {
        for (const face& f : bb.faces)
        {
            os  << 'f';
            for (const label fpi : f)
            {
                os  << ' ' << fpi + vertI + 1;
            }
            os  << nl;
        }
    }

    vertI += pts.size();
}


template<class Type>
void Foam::AABBTree<Type>::writeOBJ
(
    const bool leavesOnly,
    const bool writeLinesOnly,
    const treeBoundBox& bb,
    const label nodeI,
    const List<Pair<treeBoundBox>>& bbs,
    const List<Pair<label>>& nodes,
    label& vertI,
    Ostream& os
) const
{
    if (!leavesOnly || nodeI < 0)
    {
        writeOBJ(writeLinesOnly, bb, vertI, os);
    }

    // recurse to find leaves
    if (nodeI >= 0)
    {
        writeOBJ
        (
            leavesOnly,
            writeLinesOnly,
            bbs[nodeI].first(),
            nodes[nodeI].first(),
            bbs,
            nodes,
            vertI,
            os
        );
        writeOBJ
        (
            leavesOnly,
            writeLinesOnly,
            bbs[nodeI].second(),
            nodes[nodeI].second(),
            bbs,
            nodes,
            vertI,
            os
        );
    }
}


template<class Type>
void Foam::AABBTree<Type>::createBoxes
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
) const
{
    const vector span = bb.span();

    // Determine which direction to divide the box

    direction maxDir = 0;
    scalar maxSpan = span[maxDir];
    for (label dirI = 1; dirI < 3; ++dirI)
    {
        if (span[dirI] > maxSpan)
        {
            maxSpan = span[dirI];
            maxDir = dirI;
        }
    }


    scalar divide;

    if (equalBinSize)
    {
        // Pick up points used by this set of objects

        PackedBoolList isUsedPoint(points.size());
        DynamicList<scalar> component(points.size());

        for (const label objI : objectIDs)
        {
            const Type& obj = objects[objI];

            for (const label pointI : obj)
            {
                if (isUsedPoint.set(pointI))
                {
                    component.append(points[pointI][maxDir]);
                }
            }
        }

        // Determine the median

        Foam::sort(component);

        divide = component[component.size()/2];
    }
    else
    {
        // Geometric middle
        divide = bb.min()[maxDir] + 0.5*maxSpan;
    }


    scalar divMin = divide + tolerance_*maxSpan;
    scalar divMax = divide - tolerance_*maxSpan;


    // Assign the objects to min or max bin

    DynamicList<label> minBinObjectIDs(objectIDs.size());
    treeBoundBox minBb(boundBox::invertedBox);

    DynamicList<label> maxBinObjectIDs(objectIDs.size());
    treeBoundBox maxBb(boundBox::invertedBox);

    for (const label objI : objectIDs)
    {
        const Type& obj = objects[objI];

        bool intoMin = false;
        bool intoMax = false;

        for (const label pointI : obj)
        {
            const point& pt = points[pointI];
            if (pt[maxDir] < divMin)
            {
                intoMin = true;
            }
            if (pt[maxDir] > divMax)
            {
                intoMax = true;
            }
        }

        // Note: object is inserted into both min/max child boxes (duplicated)
        // if it crosses the bin boundaries
        if (intoMin)
        {
            minBinObjectIDs.append(objI);
            boundBox tmpBb(points, obj);
            minBb.min() = min(minBb.min(), tmpBb.min());
            minBb.max() = min(minBb.max(), tmpBb.max());
        }
        if (intoMax)
        {
            maxBinObjectIDs.append(objI);
            boundBox tmpBb(points, obj);
            maxBb.min() = min(maxBb.min(), tmpBb.min());
            maxBb.max() = min(maxBb.max(), tmpBb.max());
        }
    }

    // Inflate box in case geometry reduces to 2-D
    if (minBinObjectIDs.size())
    {
        minBb.inflate(0.01);
    }
    if (maxBinObjectIDs.size())
    {
        maxBb.inflate(0.01);
    }

    minBinObjectIDs.shrink();
    maxBinObjectIDs.shrink();


    label minI;
    if (minBinObjectIDs.size() > minLeafSize_ && level < maxLevel_)
    {
        // New leaf
        minI = nodes.size();
        nodes.append(labelPair(-1, -1));
    }
    else
    {
        // Update existing leaf
        minI = -addressing.size() - 1;
        addressing.append(minBinObjectIDs);
    }

    label maxI;
    if (maxBinObjectIDs.size() > minLeafSize_ && level < maxLevel_)
    {
        // New leaf
        maxI = nodes.size();
        nodes.append(labelPair(-1, -1));
    }
    else
    {
        // Update existing leaf
        maxI = -addressing.size() - 1;
        addressing.append(maxBinObjectIDs);
    }

    nodes(nodeI) = labelPair(minI, maxI);
    bbs(nodeI) = Pair<treeBoundBox>(minBb, maxBb);

    // Recurse
    if (minI >= 0)
    {
        createBoxes
        (
            equalBinSize,
            level + 1,
            objects,
            points,
            minBinObjectIDs,
            minBb,
            minI,
            bbs,
            nodes,
            addressing
        );
    }
    if (maxI >= 0)
    {
        createBoxes
        (
            equalBinSize,
            level + 1,
            objects,
            points,
            maxBinObjectIDs,
            maxBb,
            maxI,
            bbs,
            nodes,
            addressing
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::AABBTree<Type>::AABBTree()
:
    maxLevel_(0),
    minLeafSize_(0),
    boundBoxes_(),
    addressing_()
{}


template<class Type>
Foam::AABBTree<Type>::AABBTree
(
    const UList<Type>& objects,
    const pointField& points,
    const bool equalBinSize,
    const label maxLevel,
    const label minLeafSize
)
:
    maxLevel_(maxLevel),
    minLeafSize_(minLeafSize),
    boundBoxes_(),
    addressing_()
{
    if (objects.empty())
    {
        return;
    }


    DynamicList<Pair<treeBoundBox>> bbs(maxLevel);
    DynamicList<labelPair> nodes(maxLevel);
    DynamicList<labelList> addr(maxLevel);

    nodes.append(labelPair(-1, -1));
    treeBoundBox topBb(points);
    topBb.inflate(0.01);

    DynamicList<label> objectIDs(identity(objects.size()));

    createBoxes
    (
        equalBinSize,
        0,          // starting at top level
        objects,
        points,
        objectIDs,
        topBb,
        0,          // starting node

        bbs,
        nodes,
        addr
    );


    //{
    //    OFstream os("tree.obj");
    //    label vertI = 0;
    //    writeOBJ
    //    (
    //        true,       // leavesOnly
    //        false,      // writeLinesOnly
    //
    //        topBb,
    //        0,
    //        bbs,
    //        nodes,
    //        vertI,
    //        os
    //    );
    //}


    // transfer flattened tree to persistent storage
    DynamicList<treeBoundBox> boundBoxes(2*bbs.size());
    DynamicList<labelList> addressing(2*addr.size());

    forAll(nodes, nodeI)
    {
        if (nodes[nodeI].first() < 0)
        {
            boundBoxes.append(bbs[nodeI].first());
            addressing.append(addr[-(nodes[nodeI].first() + 1)]);
        }
        if (nodes[nodeI].second() < 0)
        {
            boundBoxes.append(bbs[nodeI].second());
            addressing.append(addr[-(nodes[nodeI].second() + 1)]);
        }
    }

    boundBoxes_.transfer(boundBoxes);
    addressing_.transfer(addressing);


    if (0)
    {
        PackedBoolList checked(objects.size());
        for (const auto& box : addressing_)
        {
            for (const auto& id : box)
            {
                checked.set(id);
            }
        }

        const label unsetSize = checked.size() - checked.count();

        if (unsetSize)
        {
            Info<< "*** Problem: IDs not set: " << unsetSize << endl;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::List<Foam::treeBoundBox>& Foam::AABBTree<Type>::boundBoxes() const
{
    return boundBoxes_;
}


template<class Type>
const Foam::List<Foam::labelList>& Foam::AABBTree<Type>::addressing() const
{
    return addressing_;
}


template<class Type>
bool Foam::AABBTree<Type>::pointInside(const point& pt) const
{
    for (const treeBoundBox& bb : boundBoxes_)
    {
        if (bb.contains(pt))
        {
            return true;
        }
    }

    return false;
}


template<class Type>
bool Foam::AABBTree<Type>::overlaps(const boundBox& bbIn) const
{
    for (const treeBoundBox& bb : boundBoxes_)
    {
        if (bb.overlaps(bbIn))
        {
            return true;
        }
    }

    return false;
}


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const AABBTree<Type>& tree)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << tree.maxLevel_ << token::SPACE
            << tree.minLeafSize_ << token::SPACE
            << tree.boundBoxes_ << token::SPACE
            << tree.addressing_ << token::SPACE;
    }
    else
    {
        os.write
        (
            reinterpret_cast<const char*>(&tree.maxLevel_),
            sizeof(tree.maxLevel_)
          + sizeof(tree.minLeafSize_)
        );
        os  << tree.boundBoxes_
            << tree.addressing_;
    }

    os.check(FUNCTION_NAME);
    return os;
}


template<class Type>
Foam::Istream& Foam::operator>>(Istream& is, AABBTree<Type>& tree)
{
    if (is.format() == IOstream::ASCII)
    {
        is  >> tree.maxLevel_
            >> tree.minLeafSize_;
    }
    else
    {
        is>> tree.maxLevel_;
        is>> tree.minLeafSize_;
//         is.beginRawRead();
//
//         readRawLabel(is, &tree.maxLevel_);
//         readRawLabel(is, &tree.minLeafSize_);
//
//         is.endRawRead();
    }

    is  >> tree.boundBoxes_
        >> tree.addressing_;

    is.check(FUNCTION_NAME);
    return is;
}


// ************************************************************************* //
