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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Get type of octant
template <class Type>
inline Foam::label Foam::treeNode<Type>::getVolType(const label octant) const
{
    return (volType_ >> 2*octant) & 0x3;
}


template <class Type>
inline const Foam::point& Foam::treeNode<Type>::midpoint() const
{
    return mid_;
}


template <class Type>
inline Foam::treeElem<Type>* const* Foam::treeNode<Type>::subNodes() const
{
    return subNodes_;
}


// octant contains pointer to treeNode(1) or treeLeaf(0)
template <class Type>
inline Foam::label Foam::treeNode<Type>::isNode(const label octant) const
{
    return subNodeTypes_ & (0x1 << octant);
}


// Get pointer to sub node
template <class Type>
inline Foam::treeNode<Type>* Foam::treeNode<Type>::getNodePtr
(
    const label octant
) const
{
#   ifdef FULLDEBUG
    if (!isNode(octant))
    {
        FatalErrorIn("treeNode::getNodePtr(const label)")
            << "not a treeNode"
            << abort(FatalError);
    }
#   endif

    return static_cast<treeNode<Type>*>(subNodes_[octant]);
}


// Get pointer to sub leaf
template <class Type>
inline Foam::treeLeaf<Type>* Foam::treeNode<Type>::getLeafPtr
(
    const label octant
) const
{
#   ifdef FULLDEBUG
    if (isNode(octant))
    {
        FatalErrorIn("treeNode::getLeafPtr(const label)")
            << "not a treeLeaf"
            << abort(FatalError);
    }
#   endif

    return static_cast<treeLeaf<Type>*>(subNodes_[octant]);
}


// ************************************************************************* //
