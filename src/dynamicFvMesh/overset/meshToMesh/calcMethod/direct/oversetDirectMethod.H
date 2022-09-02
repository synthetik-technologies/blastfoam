/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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
    Foam::oversetDirectMethod

Description
    Direct (one-to-one cell correspondence) mesh-to-mesh interpolation class

    Volume conservative.

SourceFiles
    oversetDirectMethod.C

\*---------------------------------------------------------------------------*/

#ifndef oversetDirectMethod_H
#define oversetDirectMethod_H

#include "oversetMeshToMeshMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class oversetDirectMethod Declaration
\*---------------------------------------------------------------------------*/

class oversetDirectMethod
:
    public oversetMeshToMeshMethod
{
protected:

    // Protected Member Functions

        //- Return the true if cells intersect
        virtual bool intersect
        (
            const label srcCelli,
            const label tgtCelli
        ) const;

        //- Find indices of overlapping cells in src and tgt meshes - returns
        //  true if found a matching pair
        virtual bool findInitialSeeds
        (
            const labelList& srcCellIDs,
            const boolList& mapFlag,
            const label startSeedI,
            label& srcSeedI,
            label& tgtSeedI
        ) const;

        //- Calculate the mesh-to-mesh addressing and weights
        virtual void calculateAddressing
        (
            labelListList& srcToTgtCellAddr,
            scalarListList& srcToTgtCellWght,
            labelListList& tgtToSrcCellAddr,
            scalarListList& tgtToSrcCellWght,
            const label srcSeedI,
            const label tgtSeedI,
            const labelList& srcCellIDs,
            boolList& mapFlag,
            label& startSeedI
        );

        //- Append to list of src mesh seed indices
        virtual void appendToDirectSeeds
        (
            boolList& mapFlag,
            labelList& srcTgtSeed,
            DynamicList<label>& srcSeeds,
            label& srcSeedI,
            label& tgtSeedI
        ) const;


private:

    // Private member functions

        //- No copy construct
        oversetDirectMethod(const oversetDirectMethod&) = delete;

        //- No copy assignment
        void operator=(const oversetDirectMethod&) = delete;


public:

    //- Run-time type information
    TypeName("direct");

    //- Construct from source and target meshes
    oversetDirectMethod(const polyMesh& src, const polyMesh& tgt);

    //- Destructor
    virtual ~oversetDirectMethod();


    // Member Functions

        // Evaluate

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
