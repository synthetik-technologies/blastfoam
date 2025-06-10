/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022-2025
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
    Foam::fvMeshTopoChangers::polyRefiner

SourceFiles
    fvMeshPolyRefiner.C

Description
    A fvMesh with run-time load balancing.

\*---------------------------------------------------------------------------*/

#ifndef fvMeshPolyRefiner_H
#define fvMeshPolyRefiner_H

#include "fvMeshRefiner.H"
#include "polyMeshPolyRefiner.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{

/*---------------------------------------------------------------------------*\
                Class fvMeshTopoChangers::polyRefiner Declaration
\*---------------------------------------------------------------------------*/

class polyRefiner
:
    public fvMeshRefiner
{

public:

    //- Runtime type information
    TypeName(polyMeshPolyRefiner::typeName_());


    // Constructors

        //- Construct from fvMesh and dictionary
        polyRefiner
        (
            fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        polyRefiner(const polyRefiner&) = delete;


    //- Destructor
    virtual ~polyRefiner();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const polyRefiner&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshTopoChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
