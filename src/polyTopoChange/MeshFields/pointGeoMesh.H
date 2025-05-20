/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
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
    Foam::pointGeoMesh

SourceFiles

Description
    Return number of points from the mesh

\*---------------------------------------------------------------------------*/

#ifndef pointGeoMesh_H
#define pointGeoMesh_H

#include "GeoMesh.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class pointGeoMesh Declaration
\*---------------------------------------------------------------------------*/

class pointGeoMesh
:
    public GeoMesh<polyMesh>
{
public:

    // Constructors

        //- Construct from fvMesh
        explicit pointGeoMesh(const polyMesh& mesh)
        :
            GeoMesh<polyMesh>(mesh)
        {}


    // Member Functions

        //- Return size
        label size() const
        {
            return size(mesh_);
        }

        //- Return size
        static label size(const Mesh& mesh)
        {
            return mesh.nPoints();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
