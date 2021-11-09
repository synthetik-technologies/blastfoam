/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "immersedBox.H"
#include "addToRunTimeSelectionTable.H"
#include "immersedBoundaryObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedBox, 0);
    addToRunTimeSelectionTable(immersedShape, immersedBox, threeD);
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::immersedBox::immersedBox
(
    const polyMesh& pMesh,
    const immersedBoundaryObject& ibo,
    const dictionary& dict
)
:
    immersedShape(pMesh, ibo, dict),
    L_(dict.lookup<vector>("L"))
{
    read(dict);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedBox::~immersedBox()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::autoPtr<Foam::standAlonePatch> Foam::immersedBox::createPatch() const
{
    pointField points;
    List<face> faces;
    makeBox(L_, points, faces);
    triSurface tri(triFaceList(faces), points);
    refine3D(tri);
    return autoPtr<standAlonePatch>
    (
        new standAlonePatch(tri.faces(), tri.points())
    );
}


Foam::labelList
Foam::immersedBox::calcInside(const pointField& points) const
{
    labelList insidePoints(points.size(), -1);
    pointField validPoints(points.size());
    label pi = 0;
    forAll(points, i)
    {
        if (bb_.contains(points[i]))
        {
            validPoints[pi] = points[i];
            insidePoints[pi++] = i;
        }
    }
    insidePoints.resize(pi);
    validPoints.resize(pi);
    if (!pi)
    {
        return insidePoints;
    }
    pi = 0;

    boundBox bb(-L_/2.0, L_/2.0);
    validPoints = object_.inverseTransform(validPoints);

    forAll(validPoints, i)
    {
        if (bb.contains(validPoints[i]))
        {
            insidePoints[pi++] = insidePoints[i];
        }
    }
    insidePoints.resize(pi);
    return insidePoints;
}


bool Foam::immersedBox::inside(const point& pt) const
{
    if (bb_.contains(pt))
    {
        boundBox bb(-L_/2.0, L_/2.0);
        if (bb.contains(object_.inverseTransform(pt)))
        {
            return true;
        }
    }
    return false;
}


void Foam::immersedBox::makeBox
(
    const vector& L,
    pointField& points,
    List<face>& faces
)
{
    scalar x = L[0]/2.0;
    scalar y = L[1]/2.0;
    scalar z = L[2]/2.0;
    points = pointField
    (
        {
            vector(-x, -y, -z),
            vector(x, -y, -z),
            vector(x, y, -z),
            vector(-x, y, -z),
            vector(-x, -y, z),
            vector(x, -y, z),
            vector(x, y, z),
            vector(-x, y, z)
        }
    );
    faces = List<face>
    (
        {
            triFace(1, 3, 0),
            triFace(2, 3, 1),
            triFace(0, 3, 7),
            triFace(4, 0, 7),
            triFace(6, 5, 4),
            triFace(7, 6, 4),
            triFace(5, 2, 1),
            triFace(5, 6, 2),
            triFace(4, 5, 0),
            triFace(0, 5, 1),
            triFace(2, 6, 7),
            triFace(3, 2, 7)
        }
    );
}


void Foam::immersedBox::write(Ostream& os) const
{
    writeEntry(os, "L", L_);
}


// ************************************************************************* //
