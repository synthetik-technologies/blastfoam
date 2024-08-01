/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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
    Foam::patchToPatches::nearest

Description
    Class to generate patchToPatch coupling geometry. Couples a face to the
    single nearest opposite face only.

SourceFiles
    nearest.C

\*---------------------------------------------------------------------------*/

#ifndef nearestPatchToPatchMapping_H
#define nearestPatchToPatchMapping_H

#include "nearbyPatchToPatchMapping.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace patchToPatchMappings
{

/*---------------------------------------------------------------------------*\
                Class patchToPatchMappings::nearest Declaration
\*---------------------------------------------------------------------------*/

class nearest
:
    public nearby
{
    // Private data

        //- For each source face, the distance to its coupled target face
        List<scalar> srcFaceDistances_;

        //- For each target face, the distance to its coupled source face
        List<scalar> tgtFaceDistances_;

        //- For each source point, the distance to its coupled target point
        List<scalar> srcPointDistances_;

        //- For each target point, the distance to its coupled source point
        List<scalar> tgtPointDistances_;



    // Private Member Functions

        //- Find nearest points
        void findNearestPoints
        (
            const primitivePatch& patch,
            const primitivePatch& otherPatch,
            const List<DynamicList<label>>& faceAddr,
            List<DynamicList<label>>& pointAddr,
            List<scalar>& distances
        );

        //- Intersect two faces
        virtual bool intersectFaces
        (
            const primitivePatch& srcPatch,
            const pointField& srcPts0,
            const primitivePatch& tgtPatch,
            const pointField& tgtPts0,
            const vectorField& pointNormals,
            const vectorField& pointNormals0,
            const label srcFacei,
            const label tgtFacei
        );

        //- Initialisation
        virtual void initialise
        (
            const primitivePatch& srcPatch,
            const pointField& srcPts0,
            const primitivePatch& tgtPatch,
            const pointField& tgtPts0,
            const vectorField& pointNormals,
            const vectorField& pointNormals0
        );

        //- Finalise the intersection locally. Trims the local target patch
        //  so that only parts that actually intersect the source remain.
        //  Returns new-to-old map for trimming target-associated data.
        virtual labelList finaliseLocalPoints
        (
            const primitivePatch& srcPatch,
            const pointField& srcPts0,
            const primitivePatch& tgtPatch,
            const pointField& tgtPts0,
            const vectorField& pointNormals,
            const vectorField& pointNormals0
        );
        virtual labelList finaliseLocalFaces
        (
            const primitivePatch& srcPatch,
            const pointField& srcPts0,
            const primitivePatch& tgtPatch,
            const pointField& tgtPts0,
            const vectorField& pointNormals,
            const vectorField& pointNormals0
        );

        //- Send data that resulted from an intersection between the source
        //  patch and a distributed source-local-target patch back to the
        //  original target processes
        virtual void rDistributeTgt
        (
            const primitivePatch& tgtPatch,
            const pointField& tgtPts0
        );


protected:

        //- Finalise the intersection
        virtual label finalisePoints
        (
            const primitivePatch& srcPatch,
            const pointField& srcPts0,
            const primitivePatch& tgtPatch,
            const pointField& tgtPts0,
            const vectorField& pointNormals,
            const vectorField& pointNormals0,
            const transformer& BtoA
        );
        virtual label finaliseFaces
        (
            const primitivePatch& srcPatch,
            const pointField& srcPts0,
            const primitivePatch& tgtPatch,
            const pointField& tgtPts0,
            const vectorField& pointNormals,
            const vectorField& pointNormals0,
            const transformer& BtoA
        );


public:

    //- Runtime type information
    TypeName("nearest");


    // Constructors

        //- Construct from components
        nearest
        (
            const primitivePatch& srcPatch,
            const primitivePatch& tgtPatch,
            const dictionary& dict,
            const bool needPoints,
            const bool reverse
        );


    // Destructor

        virtual ~nearest()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchToPatchMappings
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
