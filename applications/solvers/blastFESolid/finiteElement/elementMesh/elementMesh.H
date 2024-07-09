/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::elementMesh

Description
    Mesh data needed to do the Finite Volume discretisation.

Class
    Foam::isElementMesh

Description
    Supports static assertion that a template argument is of type elementMesh.

\*---------------------------------------------------------------------------*/

#ifndef elementMesh_H
#define elementMesh_H

#include "GeoMesh.H"
#include "feMesh.H"
#include "primitiveMesh.H"
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class elementMesh Declaration
\*---------------------------------------------------------------------------*/

class elementMesh
:
    public GeoMesh<fvMesh>
{

public:

    // Constructors

        //- Construct from fvMesh
        explicit elementMesh(const fvMesh& mesh)
        :
            GeoMesh<fvMesh>(mesh)
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
            return mesh.nElements();
        }

        //- Return cell centres
        const elementVectorField& C()
        {
            return mesh_.nodes();
        }
};


/*---------------------------------------------------------------------------*\
                        Class isElementMesh Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class isElementMesh
:
    public std::false_type
{};


template<>
class isElementMesh<elementMesh>
:
    public std::true_type
{};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
