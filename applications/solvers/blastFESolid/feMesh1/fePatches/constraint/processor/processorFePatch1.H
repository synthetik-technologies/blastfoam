/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::processorFePatch1

Description
    Processor patch boundary needs to be such that the ordering of
    points in the patch is the same on both sides.

    Looking at the creation of the faces on both sides of the processor
    patch they need to be identical on both sides with the normals pointing
    in opposite directions.  This is achieved by calling the reverseFace
    function in the decomposition.  It is therefore possible to re-create
    the ordering of patch points on the slave side by reversing all the
    patch faces of the owner.

SourceFiles
    processorFePatch1.C

\*---------------------------------------------------------------------------*/

#ifndef processorFePatch1_H
#define processorFePatch1_H

#include "fePatch1.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class processorFePatch Declaration
\*---------------------------------------------------------------------------*/

class processorFePatch1
:
    public fePatch1
{
    // Private Data

        const processorPolyPatch& procPolyPatch_;

public:

    //- Runtime type information
    TypeName(processorPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        processorFePatch1
        (
            const polyPatch& patch,
            const feBoundaryMesh1& bm
        );


    //- Destructor
    virtual ~processorFePatch1();


    // Member Functions

        //- Return size
        virtual label size() const
        {
            return 0;
        }

        //- Return message tag to use for communication
        virtual int tag() const
        {
            return procPolyPatch_.tag();
        }

        //- Return communicator used for communication
        virtual label comm() const
        {
            return procPolyPatch_.comm();
        }

        //- Return processor number
        int myProcNo() const
        {
            return procPolyPatch_.myProcNo();
        }

        //- Return neighbour processor number
        int neighbProcNo() const
        {
            return procPolyPatch_.neighbProcNo();
        }

        //- Is this a master patch
        bool isMaster() const
        {
            return myProcNo() < neighbProcNo();
        }

        //- Is this a slave patch
        bool isSlave() const
        {
            return !isMaster();
        }

        //- Return the underlying processorPolyPatch
        const processorPolyPatch& procPolyPatch() const
        {
            return procPolyPatch_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
