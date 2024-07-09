/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::internalFePatch1

Description
    Constraint patch to hold internal faces exposed by sub-setting.

SourceFiles
    internalFePatch1.C


\*---------------------------------------------------------------------------*/

#ifndef internalFePatch1_H
#define internalFePatch1_H

#include "fePatch1.H"
#include "internalPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class internalFePatch Declaration
\*---------------------------------------------------------------------------*/

class internalFePatch1
:
    public fePatch1
{

public:

    //- Runtime type information
    TypeName(internalPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        internalFePatch1
        (
            const polyPatch& patch,
            const feBoundaryMesh1& bm
        )
        :
            fePatch1(patch, bm)
        {}

    //- Return size
    virtual label size() const
    {
        return 0;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
