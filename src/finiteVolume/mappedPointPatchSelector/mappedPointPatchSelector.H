/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::mappedPointPatchSelector

Description
    Class used to return mapping functions using either a mappedPolyPatch
    or a mappedMovingPolyPatch


Note

SourceFiles
    mappedPointPatchSelector.C

\*---------------------------------------------------------------------------*/

#ifndef mappedPointPatchSelector_H
#define mappedPointPatchSelector_H

#include "pointPatch.H"
#include "mappedMovingPointPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class mappedPointPatchSelector Declaration
\*---------------------------------------------------------------------------*/

class mappedPointPatchSelector
{
protected:

    // Protected data

        //- Pointer to mappedMovingPatchBase
        const mappedMovingPointPatchBase* mappedMovingPatchPtr_;


public:

    //- Runtime type information
    TypeName("mappedPointPatchSelector");


    // Constructors

        //- Construct from patch
        mappedPointPatchSelector(const pointPatch&);


    //- Destructor
    virtual ~mappedPointPatchSelector();


    // Member Functions

        // Access

            //- Region to sample
            inline const word& sampleRegion() const;

            //- Patch (only if NEARESTPATCHFACE)
            inline const word& samplePatch() const;

            //- PatchGroup (only if NEARESTPATCHFACE)
            inline const word& coupleGroup() const;

            //- Return size of mappedMoving mesh/patch/boundary
            inline label sampleSize() const;

            //- Return reference to the parallel distribution map
            inline const mapDistribute& map() const;

            //- Get the region mesh
            const polyMesh& sampleMesh() const;

            //- Get the patch on the region
            const polyPatch& samplePolyPatch() const;


        // Helpers

            //- Get the sample points
            tmp<pointField> samplePoints() const;


        // Distribute

            //- Wrapper around map/interpolate data distribution
            template<class Type>
            void distribute(List<Type>& lst) const;

            //- Wrapper around map/interpolate data distribution with operation
            template<class Type, class CombineOp>
            void distribute(List<Type>& lst, const CombineOp& cop) const;

            //- Wrapper around map/interpolate data distribution
            template<class Type>
            void reverseDistribute(List<Type>& lst) const;

            //- Wrapper around map/interpolate data distribution with operation
            template<class Type, class CombineOp>
            void reverseDistribute(List<Type>& lst, const CombineOp& cop) const;


    // Static functions

        // Function to clear patches if an update has occurred
        static void clearMappedPatches(fvMesh& mesh);

        // Set the displacement field given a field
        static bool setMappedPatchDisplacement
        (
            fvMesh& ownMesh,
            const fvMesh& neiMesh
        );

        //- Check if a patch is of a mapped type
        static bool isAMappedType(const pointPatch& patch);

        //- Replace mappedFvPatches with  mappedMovingPatches
        static void replaceMappedWithMappedMoving(fvMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mappedPointPatchSelectorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mappedPointPatchSelectorTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
