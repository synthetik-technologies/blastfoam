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
    Foam::patchToPatches::nearby

Description
    Class to generate patchToPatch coupling geometry. Couples a face to the
    single nearby opposite face only.

SourceFiles
    nearby.C

\*---------------------------------------------------------------------------*/

#ifndef nearbyPatchToPatchMapping_H
#define nearbyPatchToPatchMapping_H

#include "patchToPatchMapping.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace patchToPatchMappings
{

/*---------------------------------------------------------------------------*\
                Class patchToPatchMappings::nearby Declaration
\*---------------------------------------------------------------------------*/

class nearby
:
    public patchToPatchMapping
{
protected:
    // Protected data

        //- Minum bounding box dimension
        scalar minBbDim_;

        //- For each source face, the bounding sphere
        List<Tuple2<point, scalar>> srcSpheres_;

        //- For each target face, the bounding sphere
        List<Tuple2<point, scalar>> tgtSpheres_;


    // Private Member Functions

        //- Get the bound box for a source face
        virtual treeBoundBox makeBb
        (
            const face& f,
            const pointField& pts,
            const vectorField& pns
        ) const;

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


public:

    //- Runtime type information
    TypeName("nearby");

    // Constructors

        //- Construct from components
        nearby
        (
            const primitivePatch& srcPatch,
            const primitivePatch& tgtPatch,
            const dictionary& dict,
            const bool needPoints,
            const bool reverse
        );

        //- Construct from components
        nearby
        (
            const primitivePatch& srcPatch,
            const primitivePatch& tgtPatch,
            const dictionary& dict,
            const bool needPoints,
            const bool reverse,
            const scalar minBbDim
        );


    // Destructor

        virtual ~nearby();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchToPatchMappings
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
