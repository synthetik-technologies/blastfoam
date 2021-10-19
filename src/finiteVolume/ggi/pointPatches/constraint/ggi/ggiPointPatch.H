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
    ggiPointPatch

Description
    Dummy ggi patch for post-processing.  No functionality built in!

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    ggiPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef ggiPointPatch_H
#define ggiPointPatch_H

#include "coupledFacePointPatch.H"
#include "ggiPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ggiPointPatch Declaration
\*---------------------------------------------------------------------------*/

class ggiPointPatch
:
    public coupledFacePointPatch
{
    // Private data

        //- Local reference cast into the cyclic patch
        const ggiPolyPatch& ggiPolyPatch_;


    // Private Member Functions

        //- Disallow default construct as copy
        ggiPointPatch(const ggiPointPatch&);

        //- Disallow default assignment
        void operator=(const ggiPointPatch&);


protected:

    // Protected Member Functions

        //- Initialise the calculation of the patch geometry
        virtual void initCalcGeometry(PstreamBuffers&)
        {}

        //- Calculate mesh points
        virtual void calcGeometry(PstreamBuffers&)
        {}

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&)
        {}

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&)
        {}

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&);

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&);


public:

    typedef pointBoundaryMesh BoundaryMesh;


    //- Runtime type information
    TypeName(ggiPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        ggiPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        );


    //- Destructor
    virtual ~ggiPointPatch();


    // Member Functions

        // Access

            //- Return true because this patch is coupled
            virtual bool coupled() const;


        // Communications support

            //- Return communicator used for communication
            virtual label comm() const
            {
                return ggiPolyPatch_.comm();
            }

            //- Return message tag to use for communication
            virtual int tag() const
            {
                return ggiPolyPatch_.tag();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
