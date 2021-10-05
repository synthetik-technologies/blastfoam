/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::staticBlastFvMesh

Description
    Foam::staticBlastFvMesh

SourceFiles
    staticBlastFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef staticBlastFvMesh_H
#define staticBlastFvMesh_H

#include "dynamicBlastFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class staticBlastFvMesh Declaration
\*---------------------------------------------------------------------------*/

class staticBlastFvMesh
:
    public dynamicBlastFvMesh
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        staticBlastFvMesh(const staticBlastFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const staticBlastFvMesh&);


public:

    //- Runtime type information
    TypeName("staticFvMesh");


    // Constructors

        //- Construct from IOobject
        staticBlastFvMesh(const IOobject& io);


    //- Destructor
    ~staticBlastFvMesh();


    // Member Functions

        //- Dummy update function which does not change the mesh
        virtual bool update();

        //- Refine the mesh
        virtual bool refine(const bool correctError = true);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
