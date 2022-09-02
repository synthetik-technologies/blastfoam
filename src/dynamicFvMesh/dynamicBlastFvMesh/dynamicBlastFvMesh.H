
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
    Foam::dynamicBlastFvMesh

Description
    Abstract base class for geometry and/or topology changing fvMesh.

SourceFiles
    dynamicBlastFvMesh.C
    dynamicBlastFvMeshNew.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicBlastFvMesh_H
#define dynamicBlastFvMesh_H

#include "fvMesh.H"
#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dynamicBlastFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicBlastFvMesh
:
    virtual public dynamicFvMesh
{
    // Private data

        //- The dynamicMeshDict
        IOdictionary dynamicMeshDict_;

        //- Return the dynamicMeshDict IOobject
        static IOobject dynamicMeshDictIOobject(const IOobject& io);

public:

    //- Runtime type information
    TypeName("dynamicBlastFvMesh");

    // Constructors

        //- Construct from objectRegistry, and read/write options
        explicit dynamicBlastFvMesh(const IOobject& io);

        //- Construct from components without boundary.
        //  Boundary is added using addFvPatches() member function
        dynamicBlastFvMesh
        (
            const IOobject& io,
            pointField&& points,
            faceList&& faces,
            labelList&& allOwner,
            labelList&& allNeighbour,
            const bool syncPar = true
        );

        //- Construct without boundary from cells rather than owner/neighbour.
        //  Boundary is added using addPatches() member function
        dynamicBlastFvMesh
        (
            const IOobject& io,
            pointField&& points,
            faceList&& faces,
            cellList&& cells,
            const bool syncPar = true
        );

        //- Disallow default bitwise copy construction
        dynamicBlastFvMesh(const dynamicBlastFvMesh&) = delete;


    //- Destructor
    virtual ~dynamicBlastFvMesh();


    // Member Functions

        //- Return the dynamicMeshDict
        const dictionary& dynamicMeshDict() const
        {
            return dynamicMeshDict_;
        }

        //- Refine the mesh
        virtual bool refine() = 0;
};

//- Refine is type is a dynamicBlastFvMesh
bool refineMesh(dynamicFvMesh&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
