/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    directPatchToPatchMapping

Description
    patchToPatchMapping wrapper for direct mapping i.e. conformal
    meshes

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    directPatchToPatchMapping.C

\*---------------------------------------------------------------------------*/

#ifndef directPatchToPatchMapping_H
#define directPatchToPatchMapping_H

#include "nearestPatchToPatchMapping.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace patchToPatchMappings
{

/*---------------------------------------------------------------------------*\
                        Class directPatchToPatchMapping Declaration
\*---------------------------------------------------------------------------*/

class direct
:
    public nearest
{
    // Private member functions

        // Make sure every face references exactly one face
        static bool forwardCheck
        (
            const primitivePatch& patch,
            const List<DynamicList<label>>& localOtherFaces,
            const bool isSrc
        );

        // Make sure every face is referenced by exactly one face
        static bool reverseCheck
        (
            const primitivePatch& patch,
            const List<DynamicList<label>>& otherLocalFaces,
            const autoPtr<mapDistribute>& mapPtr,
            const bool isSrc
        );

        //- Finalise the intersection
        virtual label finalise
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
    TypeName("direct");


    // Constructors

        //- Construct from components
        direct
        (
            const primitivePatch& srcPatch,
            const primitivePatch& tgtPatch,
            const dictionary& dict,
            const bool reverse
        );


    // Destructor

        virtual ~direct()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchToPatchMappings
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
