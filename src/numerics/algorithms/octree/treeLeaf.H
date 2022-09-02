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
    Foam::treeLeaf

Description
    An octree treeLeaf.

SourceFiles
    treeLeaf.C
    octreeDataPointTreaLeaf.H       (specialization for points)
    octreeDataPointTreeLeaf.C
    octreeDataTriSurfaceTreeLeaf.H  (specialization for triSurface)
    octreeDataTriSurfaceTreeLeaf.C

\*---------------------------------------------------------------------------*/

#ifndef treeLeaf_H
#define treeLeaf_H

#include "labelList.H"
#include "treeElem.H"
#include "boolList.H"
#include "linePointRef.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class treeBoundBox;
class Ostream;

template<class Type> class octree;

// Forward declaration of friend functions and operators
template<class Type> class treeLeaf;

template<class Type> Istream& operator>>(Istream&, treeLeaf<Type>&);
template<class Type> Ostream& operator<<(Ostream&, const treeLeaf<Type>&);


/*---------------------------------------------------------------------------*\
                        Class treeLeafName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(treeLeaf);


/*---------------------------------------------------------------------------*\
                           Class treeLeaf Declaration
\*---------------------------------------------------------------------------*/

template <class Type>
class treeLeaf
:
    public treeElem<Type>,
    public treeLeafName
{
    // Private data

        // Keeps real size (at construction time indices_ might be untrimmed)
        label size_;

        // Indices of 'things' whose bb overlaps leaf bb.
        labelList indices_;


    // Private Member Functions

        static void space(Ostream&, const label);

        //- Disallow default bitwise copy construct
        treeLeaf(const treeLeaf&);

        //- Disallow default bitwise assignment
        void operator=(const treeLeaf&);


public:

    // Constructors

        //- Construct with size
        treeLeaf(const treeBoundBox& bb, const label size);

        //- Construct from list
        treeLeaf(const treeBoundBox& bb, const labelList& indices);

        //- Construct from Istream
        treeLeaf(Istream&);


    // Destructor

        ~treeLeaf();


    // Member Functions

        // Access

            label size() const
            {
                return size_;
            }

            const labelList& indices() const
            {
                return indices_;
            }

        // Edit

            void insert(const label index)
            {
                if (size_ >= indices_.size())
                {
                    FatalErrorIn
                    (
                        "treeLeaf<Type>::insert(index)"
                    )
                        << "overflow"
                        << "  size_ :" << size_
                        << "  size():" << indices_.size()
                        << abort(FatalError);
                }
                indices_[size_++] = index;
            }

            void trim()
            {
                if (size_ == 0)
                {
                    FatalErrorIn
                    (
                        "treeLeaf<Type>::trim()"
                    )
                        << "Trying to trim empty leaf: " << endl
                        << "  size_ :" << size_
                        << "  size():" << indices_.size()
                        << abort(FatalError);
                }
                indices_.setSize(size_);
            }

            //- Take indices at refineLevel and distribute them to lower levels
            treeLeaf<Type>* redistribute
            (
                const label,
                octree<Type>&,
                const Type&
            );

            label setSubNodeType
            (
                const label level,
                octree<Type>& top,
                const Type& shapes
            ) const;

        // Search

            //- Get type of sample
            label getSampleType
            (
                const label level,
                const octree<Type>& top,
                const Type& shapes,
                const point& sample
            ) const;

            //- Find index of shape containing sample
            label find
            (
                const Type& shapes,
                const point& sample
            ) const;

            //- Find tightest fitting bounding box in leaf
            bool findTightest
            (
                const Type& shapes,
                const point& sample,
                treeBoundBox& tightest
            ) const;

            //- Find nearest point.
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

        // Write

            //- Debug: print a leaf
            void printLeaf(Ostream&, const label) const;

            //- Debug: Write bb in OBJ format
            void writeOBJ
            (
                Ostream& os,
                const label level,
                label& vertNo
            ) const;

            //- debug:
            label countLeaf(Ostream&, const label) const;


    // IOstream Operators

        friend Istream& operator>> <Type>(Istream&, treeLeaf<Type>&);
        friend Ostream& operator<< <Type>(Ostream&, const treeLeaf<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "treeLeaf.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "octreeDataPointTreeLeaf.H"
//#include "octreeDataTriSurfaceTreeLeaf.H"  HJ, header moved to triSurfaceMesh.H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
