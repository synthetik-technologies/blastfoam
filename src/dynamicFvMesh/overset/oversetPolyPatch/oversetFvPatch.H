/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Foam::oversetFvPatch

Description
    Patch for indicating interpolated boundaries (in overset meshes).

SourceFiles
    oversetFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef oversetFvPatch_H
#define oversetFvPatch_H

#include "fvPatch.H"
#include "oversetPolyPatch.H"
#include "fvBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class oversetFvPatch Declaration
\*---------------------------------------------------------------------------*/

class oversetFvPatch
:
    public fvPatch
{
    // Private data

        const oversetPolyPatch& oversetPolyPatch_;


public:

    //- Runtime type information
    TypeName(oversetPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        oversetFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            fvPatch(patch, bm),
            oversetPolyPatch_(refCast<const oversetPolyPatch>(patch))
        {}


    // Member functions

        // Access

            //- Return local reference cast into the overset patch
            const oversetPolyPatch& oversetPatch() const
            {
                return oversetPolyPatch_;
            }


        // Access

            //- Return faceCell addressing
            virtual const labelUList& faceCells() const
            {
                return oversetPolyPatch_.faceCells();
            }

            //- Am I the master interface
            virtual bool master() const
            {
                return oversetPolyPatch_.master();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
