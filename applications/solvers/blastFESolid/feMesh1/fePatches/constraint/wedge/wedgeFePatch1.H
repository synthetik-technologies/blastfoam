/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2025
     \\/     M anipulation  | Synthetik Applied Technologies
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
    Foam::wedgeFePatch1

Description
    Wedge front and back plane patch.

SourceFiles
    wedgeFePatch1.C

\*---------------------------------------------------------------------------*/

#ifndef wedgeFePatch1_H
#define wedgeFePatch1_H

#include "fePatch1.H"
#include "wedgePolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class wedgePointPatch Declaration
\*---------------------------------------------------------------------------*/

class wedgeFePatch1
:
    public fePatch1
{
    // Private Data

        //- Local reference cast into the symmetryPlane patch
        const wedgePolyPatch& wedgePolyPatch_;


public:

    //- Runtime type information
    TypeName(wedgePolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        wedgeFePatch1
        (
            const polyPatch& patch,
            const feBoundaryMesh1& bm
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
