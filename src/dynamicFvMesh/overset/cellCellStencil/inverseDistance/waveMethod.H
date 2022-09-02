/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::waveMethod

Description
    Direct (one-to-one cell correspondence) mesh-to-mesh interpolation class

SourceFiles
    waveMethod.C

\*---------------------------------------------------------------------------*/

#ifndef waveMethod_H
#define waveMethod_H

#include "oversetMeshToMeshMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class waveMethod Declaration
\*---------------------------------------------------------------------------*/

class waveMethod
:
    public oversetMeshToMeshMethod
{
    // Private member functions

        //- No copy construct
        waveMethod(const waveMethod&) = delete;

        //- No copy assignment
        void operator=(const waveMethod&) = delete;


public:

    //- Run-time type information
    TypeName("wave");

    //- Construct from source and target meshes
    waveMethod(const polyMesh& src, const polyMesh& tgt);

    //- Destructor
    virtual ~waveMethod();


    // Member Functions

        // Evaluate

            //- Calculate addressing
            static void calculate
            (
                const polyMesh& src,
                const polyMesh& tgt,
                labelList& srcToTgtAddr
            );

            //- Calculate addressing and weights and optionally offset vectors
            virtual void calculate
            (
                labelListList& srcToTgtAddr,
                scalarListList& srcToTgtWght,
                List<List<point>>& srcToTgtVec,
                labelListList& tgtToSrcAddr,
                scalarListList& tgtToSrcWght,
                List<List<point>>& tgtToSrcVec
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
