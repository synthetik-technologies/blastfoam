/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "refineCrackerFvMesh.H"
#include "crackerFvMesh.H"
#include "simpleCrackerFvMesh.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef RefineCrackerFvMesh<crackerFvMesh> refineCrackerFvMesh;
defineTemplateTypeNameAndDebug
(
    refineCrackerFvMesh,
    0
);
addNamedToRunTimeSelectionTable
(
    topoChangerFvMesh,
    refineCrackerFvMesh,
    IOobject,
    refineCrackerFvMesh
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef RefineCrackerFvMesh<simpleCrackerFvMesh> simpleRefineCrackerFvMesh;
defineTemplateTypeNameAndDebug
(
    simpleRefineCrackerFvMesh,
    0
);
addNamedToRunTimeSelectionTable
(
    topoChangerFvMesh,
    simpleRefineCrackerFvMesh,
    IOobject,
    simpleRefineCrackerFvMesh
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
