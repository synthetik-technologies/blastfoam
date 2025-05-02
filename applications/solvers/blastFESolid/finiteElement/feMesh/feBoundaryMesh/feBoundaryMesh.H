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
    Foam::feBoundaryMesh

Description
    Foam::feBoundaryMesh

SourceFiles
    feBoundaryMesh.C

\*---------------------------------------------------------------------------*/

#ifndef feBoundaryMesh_H
#define feBoundaryMesh_H

#include "fePatchList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class feMesh;
class polyBoundaryMesh;

/*---------------------------------------------------------------------------*\
                       Class feBoundaryMesh Declaration
\*---------------------------------------------------------------------------*/

class feBoundaryMesh
:
    public fePatchList
{
    // Private Data

        //- Reference to mesh
        const feMesh& mesh_;


    // Private Member Functions

        //- Calculate the geometry for the patches (transformation tensors etc.)
        void calcGeometry();


public:

    friend class feMesh;


    // Constructors

        //- Construct from polyBoundaryMesh
        feBoundaryMesh(const feMesh&, const polyBoundaryMesh&);

        //- Disallow default bitwise copy construction
        feBoundaryMesh(const feBoundaryMesh&);


    // Member Functions

        // Access

            //- Return the mesh reference
            const feMesh& mesh() const
            {
                return mesh_;
            }

            //- Find patch index given a name
            label findPatchID(const word& patchName) const;

            //- Find patch indices given a name
            labelList findIndices(const wordRe&, const bool useGroups) const;


        //- Correct patches after moving points
        void movePoints(const pointField& p);

        //- Correct pointBoundaryMesh after topology update
        void updateMesh();

        //- Create pointBoundaryMesh from polyBoundaryMesh
        void reset(const polyBoundaryMesh&);

        //- Reorders patches. Ordering does not have to be done in
        //  ascending or descending order. Reordering has to be unique.
        //  (is shuffle). If validBoundary call needs to be parallel
        //  sync
        void shuffle(const labelUList& newToOld, const bool validBoundary);


    // Member Operators

        //- Return const and non-const reference to fvPatch by index.
        using fePatchList::operator[];

        //- Return const reference to fvPatch by name.
        const fePatch& operator[](const word&) const;

        //- Return reference to fvPatch by name.
        fePatch& operator[](const word&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const feBoundaryMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
