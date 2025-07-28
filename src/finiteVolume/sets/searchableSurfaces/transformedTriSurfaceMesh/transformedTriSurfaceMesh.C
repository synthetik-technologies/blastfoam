/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "transformedTriSurfaceMesh.H"
#include "unitConversion.H"
#include "transformer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(transformedTriSurfaceMesh, 0);
    addToRunTimeSelectionTable(searchableSurface, transformedTriSurfaceMesh, dict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::transformedTriSurfaceMesh::transformedTriSurfaceMesh
(
    const IOobject& io,
    const triSurface& s
)
:
    triSurfaceMesh(io, s)
{}


Foam::transformedTriSurfaceMesh::transformedTriSurfaceMesh(const IOobject& io)
:
    triSurfaceMesh(io)
{}


Foam::transformedTriSurfaceMesh::transformedTriSurfaceMesh
(
    const IOobject& io,
    const dictionary& dict
)
:
    triSurfaceMesh(io, dict)
{
    if (dict.isDict("transforms"))
    {
        transformer transforms;
        const dictionary& transformDict = dict.subDict("transforms");
        if (transformDict.found("rotate"))
        {
            Pair<vector> n1n2(transformDict.lookup("rotate"));

            n1n2[0] /= mag(n1n2[0]);
            n1n2[1] /= mag(n1n2[1]);

            transforms =
                transformer::rotation(rotationTensor(n1n2[0], n1n2[1]))
              & transforms;
        }
        else if (transformDict.found("Rx"))
        {
            const scalar a(transformDict.lookup<scalar>("Rx"));
            transforms =
                transformer::rotation(Rx(degToRad(a)))
              & transforms;
        }
        else if (transformDict.found("Ry"))
        {
            const scalar a(transformDict.lookup<scalar>("Ry"));
            transforms =
                transformer::rotation(Ry(degToRad(a)))
              & transforms;
        }
        else if (transformDict.found("Rz"))
        {
            const scalar a(transformDict.lookup<scalar>("Rz"));
            transforms =
                transformer::rotation(Rz(degToRad(a)))
              & transforms;
        }
        else if (transformDict.found("Ra"))
        {
            ITstream is(transformDict.lookup("Ra"));
            const vector v(is);
            const scalar a(readScalar(is));
            transforms =
                transformer::rotation(Ra(v, degToRad(a)))
              & transforms;
        }
        if (transformDict.found("scale"))
        {
            const vector v(transformDict.lookup("scale"));
            transforms =
                transformer::scaling(diagTensor(v.x(), v.y(), v.z()))
              & transforms;
        }
        if (transformDict.found("translate"))
        {
            const vector v(transformDict.lookup("translate"));
            transforms = transformer::translation(v) & transforms;
        }

        if (transforms.transformsPosition())
        {
            this->setPoints(transforms.transformPosition(triSurface::points()));

            bounds() = boundBox(triSurface::points());
        }
    }
}


Foam::transformedTriSurfaceMesh::transformedTriSurfaceMesh
(
    const IOobject& io,
    const bool isGlobal
)
:
    triSurfaceMesh(io, isGlobal)
{}


Foam::transformedTriSurfaceMesh::transformedTriSurfaceMesh
(
    const IOobject& io,
    const dictionary& dict,
    const bool isGlobal
)
:
    triSurfaceMesh(io, dict, isGlobal)
{
    if (dict.isDict("transforms"))
    {
        transformer transforms;
        const dictionary& transformDict = dict.subDict("transforms");
        if (transformDict.found("rotate"))
        {
            Pair<vector> n1n2(transformDict.lookup("rotate"));

            n1n2[0] /= mag(n1n2[0]);
            n1n2[1] /= mag(n1n2[1]);

            transforms =
                transformer::rotation(rotationTensor(n1n2[0], n1n2[1]))
              & transforms;
        }
        else if (transformDict.found("Rx"))
        {
            const scalar a(transformDict.lookup<scalar>("Rx"));
            transforms =
                transformer::rotation(Rx(degToRad(a)))
              & transforms;
        }
        else if (transformDict.found("Ry"))
        {
            const scalar a(transformDict.lookup<scalar>("Ry"));
            transforms =
                transformer::rotation(Ry(degToRad(a)))
              & transforms;
        }
        else if (transformDict.found("Rz"))
        {
            const scalar a(transformDict.lookup<scalar>("Rz"));
            transforms =
                transformer::rotation(Rz(degToRad(a)))
              & transforms;
        }
        else if (transformDict.found("Ra"))
        {
            ITstream is(transformDict.lookup("Ra"));
            const vector v(is);
            const scalar a(readScalar(is));
            transforms =
                transformer::rotation(Ra(v, degToRad(a)))
              & transforms;
        }
        if (transformDict.found("scale"))
        {
            const vector v(transformDict.lookup("scale"));
            transforms =
                transformer::scaling(diagTensor(v.x(), v.y(), v.z()))
              & transforms;
        }
        if (transformDict.found("translate"))
        {
            const vector v(transformDict.lookup("translate"));
            transforms = transformer::translation(v) & transforms;
        }

        this->setPoints(transforms.transformPosition(triSurface::points()));

        bounds() = boundBox(triSurface::points());
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::transformedTriSurfaceMesh::~transformedTriSurfaceMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
