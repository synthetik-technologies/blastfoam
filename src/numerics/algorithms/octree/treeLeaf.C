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

#include "treeLeaf.H"
#include "treeNode.H"
#include "treeBoundBox.H"
#include "octree.H"
#include "HashSet.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class Type>
void Foam::treeLeaf<Type>::space(Ostream& os, const label n)
{
    for (label i=0; i<n; i++)
    {
        os<< ' ';
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct with given size
template <class Type>
Foam::treeLeaf<Type>::treeLeaf(const treeBoundBox& bb, const label size)
:
    treeElem<Type>(bb), size_(0), indices_(size)
{}


// Construct from list
template <class Type>
Foam::treeLeaf<Type>::treeLeaf(const treeBoundBox& bb, const labelList& indices)
:
    treeElem<Type>(bb), size_(indices.size()), indices_(indices)
{
}


// Construct from Istream
template <class Type>
Foam::treeLeaf<Type>::treeLeaf(Istream& is)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class Type>
Foam::treeLeaf<Type>::~treeLeaf()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Take cells at this level and distribute them to lower levels
template <class Type>
Foam::treeLeaf<Type>* Foam::treeLeaf<Type>::redistribute
(
    const label level,
    octree<Type>& top,
    const Type& shapes
)
{
    if (debug & 1)
    {
        space(Pout, level);
        Pout<< "treeLeaf::redistribute with bb:" << this->bb() << endl;
    }

    if (size_ <= top.maxLeafRatio())
    {
        // leaf small enough
        if (debug & 1)
        {
            space(Pout, level);
            Pout<< "end of treeLeaf::redistribute : small enough" << endl;
        }
        return this;
    }
    else
    {
        // create treeNode for this level
        treeNode<Type>* treeNodePtr = new treeNode<Type>(this->bb());

        top.setNodes(top.nNodes() + 1);

        treeNodePtr->distribute
        (
            level,
            top,
            shapes,
            indices_
        );

        if (debug & 1)
        {
            space(Pout, level);
            Pout<< "end of treeLeaf::redistribute : done creating node"
                << this->bb() << endl;
        }

        // return pointer to let level above know.
        return reinterpret_cast<treeLeaf<Type>*>(treeNodePtr);
    }
}


// Set type of subnodes. Since contains elements return mixed type always.
template <class Type>
Foam::label Foam::treeLeaf<Type>::setSubNodeType
(
    const label level,
    octree<Type>& top,
    const Type& shapes
) const
{
    if (size() == 0)
    {
        FatalErrorIn
        (
            "treeLeaf<Type>::setSubNodeType(const label, octree<Type>&, "
            "const Type&)"
        )   << "empty leaf. bb:" << this->bb()
            << abort(FatalError);
    }
    return octree<Type>::MIXED;
}


template <class Type>
Foam::label Foam::treeLeaf<Type>::getSampleType
(
    const label level,
    const octree<Type>& top,
    const Type& shapes,
    const point& sample
) const
{
    return shapes.getSampleType(top, sample);
}


template <class Type>
Foam::label Foam::treeLeaf<Type>::find
(
    const Type& shapes,
    const point& sample
) const
{
    forAll(indices_, i)
    {
        if (shapes.contains(indices_[i], sample))
        {
            return indices_[i];
        }
    }

    return -1;
}


template <class Type>
bool Foam::treeLeaf<Type>::findTightest
(
    const Type& shapes,
    const point& sample,
    treeBoundBox& tightest
) const
{
    bool changed = false;

    forAll(indices_, i)
    {
        changed |= shapes.findTightest
        (
            indices_[i],
            sample,
            tightest
        );
    }

    return changed;
}


template <class Type>
bool Foam::treeLeaf<Type>::findNearest
(
    const Type& shapes,
    const point& sample,
    treeBoundBox& tightest,
    label& tightestI,
    scalar& tightestDist
) const
{
    bool changed = false;

    forAll(indices_, i)
    {
        if (shapes.overlaps(indices_[i], tightest))
        {
            if (debug & 8)
            {
                //space(Pout, level);
                Pout<< "treeLeaf<Type>::findNearest : sample:" << sample
                    << "  shape:" << indices_[i] << " overlaps:" << tightest
                    << endl;
            }
            point nearest;
            scalar thisDist = shapes.calcNearest(indices_[i], sample, nearest);

            if (thisDist < tightestDist)
            {
                // Construct new tightest Bb
                point dist(thisDist, thisDist, thisDist);

                tightest.min() = sample - dist;
                tightest.max() = sample + dist;

                // Update other return values
                tightestI = indices_[i];

                tightestDist = thisDist;

                changed = true;

                if (debug & 8)
                {
                    //space(Pout, level);
                    Pout<< "treeLeaf<Type>::findNearest : Found nearer : shape:"
                        << tightestI << "  distance:" << tightestDist
                        << " to sample:" << sample << endl;
                }
            }
        }
    }

    if (changed)
    {
        if (debug & 8)
        {
            //space(Pout, level);
            Pout<< "treeLeaf<Type>::findNearest : sample:" << sample
                << "  new nearer:" << tightestDist
                << endl;
        }
    }
    return changed;
}


template <class Type>
bool Foam::treeLeaf<Type>::findNearest
(
    const Type& shapes,
    const linePointRef& ln,
    treeBoundBox& tightest,
    label& tightestI,
    point& linePoint,   // nearest point on line
    point& shapePoint   // nearest point on shape
) const
{
    // Initial smallest distance
    scalar tightestDist = mag(linePoint - shapePoint);

    bool changed = false;

    forAll(indices_, i)
    {
        if (shapes.overlaps(indices_[i], tightest))
        {
            // Calculate nearest point on line and on shape.
            point linePt, shapePt;
            scalar thisDist = shapes.calcNearest
            (
                indices_[i],
                ln,
                linePt,
                shapePt
            );

            if (thisDist < tightestDist)
            {
                // Found nearer. Use.
                tightestDist = thisDist;
                tightestI = indices_[i];
                linePoint = linePt;
                shapePoint = shapePt;
                // Construct new tightest Bb. Nearest point can never be further
                // away than bounding box of line + margin equal to the distance
                vector span(thisDist, thisDist, thisDist);

                tightest.min() = min(ln.start(), ln.end()) - span;
                tightest.max() = max(ln.start(), ln.end()) + span;

                changed = true;
            }
        }
    }

    return changed;
}


template <class Type>
bool Foam::treeLeaf<Type>::findBox
(
    const Type& shapes,
    const treeBoundBox& box,
    labelHashSet& elements
) const
{
    bool changed = false;

    forAll(indices_, i)
    {
        if (shapes.overlaps(indices_[i], box))
        {
            elements.insert(indices_[i]);

            changed = true;
        }
    }

    return changed;
}


template <class Type>
void Foam::treeLeaf<Type>::printLeaf
(
    Ostream& os,
    const label level
) const
{
    space(os, level);

    os  << "leaf:" << this->bb()
        << "   number of entries:" << indices().size() << endl;

    space(os, level);

    os << indices() << endl;
}


// Dump cube coordinates in OBJ format
template <class Type>
void Foam::treeLeaf<Type>::writeOBJ
(
    Ostream& os,
    const label level,
    label& vertNo
) const
{
    point min = this->bb().min();
    point max = this->bb().max();

    os << "v " << min.x() << " " << min.y() << " " << min.z() << endl;
    os << "v " << max.x() << " " << min.y() << " " << min.z() << endl;
    os << "v " << max.x() << " " << max.y() << " " << min.z() << endl;
    os << "v " << min.x() << " " << max.y() << " " << min.z() << endl;

    os << "v " << min.x() << " " << min.y() << " " << max.z() << endl;
    os << "v " << max.x() << " " << min.y() << " " << max.z() << endl;
    os << "v " << max.x() << " " << max.y() << " " << max.z() << endl;
    os << "v " << min.x() << " " << max.y() << " " << max.z() << endl;

    os << "l " << vertNo   << " " << vertNo+1 << endl;
    os << "l " << vertNo+1 << " " << vertNo+2 << endl;
    os << "l " << vertNo+2 << " " << vertNo+3 << endl;
    os << "l " << vertNo+3 << " " << vertNo   << endl;

    os << "l " << vertNo+4 << " " << vertNo+5 << endl;
    os << "l " << vertNo+5 << " " << vertNo+6 << endl;
    os << "l " << vertNo+6 << " " << vertNo+7 << endl;
    os << "l " << vertNo+7 << " " << vertNo   << endl;

    os << "l " << vertNo   << " " << vertNo+4 << endl;
    os << "l " << vertNo+1 << " " << vertNo+5 << endl;
    os << "l " << vertNo+2 << " " << vertNo+6 << endl;
    os << "l " << vertNo+3 << " " << vertNo+7 << endl;

    vertNo += 8;
}


template <class Type>
Foam::label Foam::treeLeaf<Type>::countLeaf
(
    Ostream& os,
    const label level
) const
{
    label nItems = size();

    space(os, level);

    os << "leaf:" << this->bb() << " has size:" << nItems << endl;

    return nItems;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template <class Type>
Foam::Istream& Foam::operator>> (Istream& is, treeLeaf<Type>& leaf)
{
    is >> leaf.bb() >> leaf.indices_;

    // Was written trimmed
    leaf.size_ = leaf.indices_.size();
    return is;
}


template <class Type>
Foam::Ostream& Foam::operator<< (Ostream& os, const treeLeaf<Type>& leaf)
{
    os << leaf.bb();

    if (leaf.indices().size() == leaf.size())
    {
        os << leaf.indices();
    }
    else
    {
        // Storage not trimmed
        os << token::SPACE << leaf.size() << token::SPACE << token::BEGIN_LIST;
        for (label i = 0; i < leaf.size(); i++)
        {
            os << token::SPACE << leaf.indices()[i];
        }
        os << token::END_LIST;
    }
    return os;
}


// ************************************************************************* //
