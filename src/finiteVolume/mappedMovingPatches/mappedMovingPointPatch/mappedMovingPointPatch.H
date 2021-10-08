/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-08-21 Synthetik Applied Technologies: Mapping of point patches
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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
    Foam::mappedMovingPointPatch

Description
    Determines a mapping between patch points and processors they're
    on.

SourceFiles
    mappedMovingPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef mappedMovingPointPatch_H
#define mappedMovingPointPatch_H

#include "facePointPatch.H"
#include "mappedMovingPolyPatch.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pointMesh;

/*---------------------------------------------------------------------------*\
                      Class mappedMovingPointPatch Declaration
\*---------------------------------------------------------------------------*/

class mappedMovingPointPatch
:
    public facePointPatch
{

public:

    //- Runtime type information
    TypeName(mappedMovingPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        mappedMovingPointPatch
        (
            const polyPatch& pp,
            const pointBoundaryMesh& bm
        )
        :
            facePointPatch(pp, bm)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
