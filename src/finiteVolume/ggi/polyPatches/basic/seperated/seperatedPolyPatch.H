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
    Foam::seperatedPolyPatch

Description
    The seperatedPolyPatch is an abstract base class for patches that
    couple regions of the computational domain e.g. cyclic and processor-processor with a seperation
    links.

SourceFiles
    seperatedPolyPatch.C
    seperatedPolyPatchMorph.C

\*---------------------------------------------------------------------------*/

#ifndef seperatedPolyPatch_H
#define seperatedPolyPatch_H

#include "coupledPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class seperatedPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class seperatedPolyPatch
:
    public coupledPolyPatch
{
protected:

    // Protected data

        // Alowing derived classes to create their own transformation
        // HJ, 24/Oct/2008

        //- Offset (distance) vector from one side of the couple to the other
        mutable vectorField separation_;

        //- Face transformation tensor
        mutable tensorField forwardT_;

        //- Neighbour-cell transformation tensor
        mutable tensorField reverseT_;


    // Protected Member Functions

        //- Calculate the transformation tensors
        //  smallDist : matching distance per face
        //  absTol    : absolute error in normal
        void calcTransformTensors
        (
            const vectorField& Cf,
            const vectorField& Cr,
            const vectorField& nf,
            const vectorField& nr,
            const scalarField& smallDist,
            const scalar absTol = defaultMatchTol_
        ) const;


public:

    //- Runtime type information
    TypeName("seperated");


    // Constructors

        //- Construct from components
        seperatedPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct from dictionary
        seperatedPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        seperatedPolyPatch(const seperatedPolyPatch&, const polyBoundaryMesh&);

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        seperatedPolyPatch
        (
            const seperatedPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );


    //- Destructor
    virtual ~seperatedPolyPatch();


    // Member Functions

        // Access

            //- Return true if patch is seperated
            virtual bool seperated() const
            {
                return true;
            }

            //- Are the seperated planes separated
            bool separated() const
            {
                return separation_.size();
            }

            //- Return the offset (distance) vector from one side of the couple
            //  to the other
            const vectorField& separation() const
            {
                return separation_;
            }


            //- Are the cyclic planes parallel
            bool parallel() const
            {
                return forwardT_.empty();
            }

            //- Return face transformation tensor
            const tensorField& forwardT() const
            {
                return forwardT_;
            }

            //- Return neighbour-cell transformation tensor
            const tensorField& reverseT() const
            {
                return reverseT_;
            }


        //- Synchronise communications of ordering for primitivePatch
        //  Used in cases when no topological change happens locally,
        //  but is happening on other processors
        virtual void syncOrder() const = 0;

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
