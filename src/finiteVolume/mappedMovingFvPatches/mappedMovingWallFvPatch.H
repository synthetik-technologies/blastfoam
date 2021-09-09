/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
21-08-2021 Synthetik Applied Technologies: |   Added lookup of displacement field
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
    Foam::mappedMovingWallFvPatch

Description
    Foam::mappedMovingWallFvPatch

SourceFiles
    mappedMovingWallFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef mappedMovingWallFvPatch_H
#define mappedMovingWallFvPatch_H

#include "wallFvPatch.H"
#include "mappedMovingWallPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class mappedMovingWallFvPatch Declaration
\*---------------------------------------------------------------------------*/

class mappedMovingWallFvPatch
:
    public wallFvPatch
{

public:

    //- Runtime type information
    TypeName(mappedMovingWallPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        mappedMovingWallFvPatch
        (
            const polyPatch& patch,
            const fvBoundaryMesh& bm
        )
        :
            wallFvPatch(patch, bm)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
