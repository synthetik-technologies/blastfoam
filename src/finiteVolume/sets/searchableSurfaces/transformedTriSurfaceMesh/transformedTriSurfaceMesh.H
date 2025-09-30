/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::transformedTriSurfaceMesh

Description
    IOoject and searching on triSurface

    Note: when constructing from dictionary has optional parameters:
        - scale     : scaling factor.
        - tolerance : relative tolerance for doing intersections
                      (see triangle::intersection)
        - minQuality: discard triangles with low quality when getting normal

SourceFiles
    transformedTriSurfaceMesh.C

\*---------------------------------------------------------------------------*/

#ifndef transformedTriSurfaceMesh_H
#define transformedTriSurfaceMesh_H

#include "triSurfaceMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class transformedTriSurfaceMesh Declaration
\*---------------------------------------------------------------------------*/

class transformedTriSurfaceMesh
:
    public triSurfaceMesh
{
public:

    //- Runtime type information
    TypeName("transformedTriSurfaceMesh");


    // Constructors

        //- Construct from triSurface
        transformedTriSurfaceMesh(const IOobject&, const triSurface&);

        //- Construct read
        transformedTriSurfaceMesh(const IOobject& io);

        //- Construct from IO and dictionary (used by searchableSurface).
        //  Dictionary may contain a 'scale' entry (eg, 0.001: mm -> m)
        transformedTriSurfaceMesh
        (
            const IOobject& io,
            const dictionary& dict
        );


        // Special constructors for use by distributedTriSurface. File search
        // status (local/global) supplied.

            transformedTriSurfaceMesh(const IOobject& io, const bool isGlobal);

            transformedTriSurfaceMesh
            (
                const IOobject& io,
                const dictionary& dict,
                const bool isGlobal
            );

        //- Disallow default bitwise copy construction
        transformedTriSurfaceMesh(const transformedTriSurfaceMesh&) = delete;


    //- Destructor
    virtual ~transformedTriSurfaceMesh();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const transformedTriSurfaceMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
