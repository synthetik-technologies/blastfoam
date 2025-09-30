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
    Foam::patchToPatches::inverseDistance

Description
    Class to generate patchToPatch coupling geometry. Couples a face to the
    single inverseDistance opposite face only.

SourceFiles
    inverseDistance.C

\*---------------------------------------------------------------------------*/

#ifndef inverseDistancePatchToPatchMapping_H
#define inverseDistancePatchToPatchMapping_H

#include "nearbyPatchToPatchMapping.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace patchToPatchMappings
{

/*---------------------------------------------------------------------------*\
                Class patchToPatchMappings::inverseDistance Declaration
\*---------------------------------------------------------------------------*/

class inverseDistance
:
    public nearby
{
    // Private Member Functions

        //- Check if a ray hits a given face
        static bool rayHitsFace
        (
            const point& p,
            const vector& r,
            const face& f,
            const pointField& ps
        );

        //- Generate face weights
        static void generateFaceWeights
        (
            const primitivePatch& patch,
            const primitivePatch& otherPatch,
            const bool reverse,
            List<DynamicList<label>>& otherFaces,
            List<DynamicList<scalar>>& weights
        );

        //- Generate point weights
        static void generatePointWeights
        (
            const primitivePatch& patch,
            const primitivePatch& otherPatch,
            const List<DynamicList<label>>& faceAddr,
            List<DynamicList<label>>& pointAddr,
            List<DynamicList<scalar>>& weights
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

        //- Generate point weights
        void generatePointWeights
        (
            const primitivePatch& srcPatch,
            const primitivePatch& tgtPatch
        );

        //- Generate face weights
        void generateFaceWeights
        (
            const primitivePatch& srcPatch,
            const primitivePatch& tgtPatch
        );

        //- Send data that resulted from an intersection between the source
        //  patch and a distributed source-local-target patch back to the
        //  original target processes
        virtual void rDistributeTgt
        (
            const primitivePatch& tgtPatch,
            const pointField& tgtPts0
        );

        //- Finalise the intersection
        virtual label finalisePoints
        (
            const primitivePatch& srcPatch,
            const pointField& srcPts0,
            const primitivePatch& tgtPatch,
            const pointField& tgtPts0,
            const vectorField& pointNormals,
            const vectorField& pointNormals0,
            const transformer& tgtToSrc
        );
        virtual label finaliseFaces
        (
            const primitivePatch& srcPatch,
            const pointField& srcPts0,
            const primitivePatch& tgtPatch,
            const pointField& tgtPts0,
            const vectorField& pointNormals,
            const vectorField& pointNormals0,
            const transformer& tgtToSrc
        );


public:

    //- Runtime type information
    TypeName("inverseDistance");


    // Constructors

        //- Construct from components
        inverseDistance
        (
            const primitivePatch& srcPatch,
            const primitivePatch& tgtPatch,
            const dictionary& dict,
            const bool needPoints,
            const bool reverse
        );


    // Destructor

        virtual ~inverseDistance()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchToPatchMappings
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
