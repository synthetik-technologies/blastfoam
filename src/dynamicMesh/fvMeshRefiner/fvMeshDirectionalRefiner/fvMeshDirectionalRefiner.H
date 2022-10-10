/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    Foam::fvMeshDirectionalRefiner

SourceFiles
    fvMeshDirectionalRefiner.C

Description
    Mesh refiner that refines in one direction (x, y, or z) using the
    undoableMeshCutter mesh refiner. This is a simplified implementation of
    the refinemMesh utility. Currently there is no ability to unrefine the
    mesh

\*---------------------------------------------------------------------------*/

#ifndef fvMeshDirectionalRefiner_H
#define fvMeshDirectionalRefiner_H

#include "fvMeshRefiner.H"
#include "polyMeshDirectionalRefiner.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class fvMeshDirectionalRefiner Declaration
\*---------------------------------------------------------------------------*/

class fvMeshDirectionalRefiner
:
    public fvMeshRefiner
{

public:

    //- Runtime type information
    TypeName(polyMeshDirectionalRefiner::typeName_());


    // Constructors

        //- Construct from fvMesh
        explicit fvMeshDirectionalRefiner(fvMesh& mesh);

        //- Construct from fvMesh and dictionary
        fvMeshDirectionalRefiner
        (
            fvMesh& mesh,
            const dictionary& dict,
            const bool force = false,
            const bool read = true
        );

        //- Disallow default bitwise copy construction
        fvMeshDirectionalRefiner(const fvMeshDirectionalRefiner&) = delete;


    //- Destructor
    virtual ~fvMeshDirectionalRefiner();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvMeshDirectionalRefiner&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
