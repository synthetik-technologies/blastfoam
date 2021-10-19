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
    Foam::octreeDataFace

Description
    Holds data for octree to work on mesh faces.

    For example, calculate (in calcNearest) the correct intersection point
    with a face.

SourceFiles
    octreeDataFace.C

\*---------------------------------------------------------------------------*/

#ifndef octreeDataFace_H
#define octreeDataFace_H

#include "treeBoundBoxList.H"
#include "faceList.H"
#include "point.H"
#include "className.H"
#include "linePointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class primitiveMesh;
template<class Type> class octree;
class polyPatch;

/*---------------------------------------------------------------------------*\
                       Class octreeDataFace Declaration
\*---------------------------------------------------------------------------*/

class octreeDataFace
{
    // Static data

        //- tolerance on linear dimensions
        static scalar tol;


    // Private data

        //- the mesh
        const primitiveMesh& mesh_;

        //- labels (in mesh indexing) of faces
        labelList meshFaces_;

        //- bbs for all above faces
        treeBoundBoxList allBb_;


    // Private Member Functions

        //- Set allBb to tight fitting bounding box
        void calcBb();

public:

    // Declare name of the class and its debug switch
    ClassName("octreeDataFace");

    // Constructors

        //- Construct from selected mesh faces.
        octreeDataFace
        (
            const primitiveMesh&,
            const labelList& meshFaces,
            const treeBoundBoxList&
        );

        //- Construct from selected mesh faces. Tight fitting bounding boxes
        //  generated internally.
        octreeDataFace
        (
            const primitiveMesh&,
            const labelList& meshFaces
        );

        //- Construct from selected mesh faces.
        octreeDataFace
        (
            const primitiveMesh&,
            const UList<const labelList*>&,
            const UList<const treeBoundBoxList*>&
        );

        //- Construct from selected mesh faces.
        //  Tight-fitting bounding boxes generated internally.
        octreeDataFace(const primitiveMesh&, const UList<const labelList*>&);

        //- Construct from all faces in patch.
        //  Tight-fitting bounding boxes generated internally.
        octreeDataFace(const polyPatch&);

        //- Construct from all boundary faces.
        //  Tight-fitting bounding boxes generated internally.
        octreeDataFace(const primitiveMesh&);

        //- Construct as copy
        octreeDataFace(const octreeDataFace&);


    // Destructor

        ~octreeDataFace();


    // Member Functions

        // Access

            const primitiveMesh& mesh() const
            {
                return mesh_;
            }

            const labelList& meshFaces() const
            {
                return meshFaces_;
            }

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
                const octree<octreeDataFace>&,
                const point&
            ) const;

            //- Does (bb of) shape at index overlap bb
            bool overlaps
            (
                const label index,
                const treeBoundBox& sampleBb
            ) const;

            //- Does shape at index contain sample
            bool contains(const label index, const point& sample) const;

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
