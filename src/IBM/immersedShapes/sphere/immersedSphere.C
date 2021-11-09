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

#include "immersedSphere.H"
#include "addToRunTimeSelectionTable.H"
#include "immersedBoundaryObject.H"
#include "constants.H"
#include "immersedBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedSphere, 0);
    addToRunTimeSelectionTable(immersedShape, immersedSphere, threeD);
}

// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::immersedSphere::immersedSphere
(
    const polyMesh& pMesh,
    const immersedBoundaryObject& ibo,
    const dictionary& dict
)
:
    immersedShape(pMesh, ibo, dict),
    radius_(dict.lookup<scalar>("radius"))
{
    read(dict);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedSphere::~immersedSphere()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::standAlonePatch>
Foam::immersedSphere::createPatch() const
{
    pointField points;
    List<face> faces;
    immersedBox::makeBox(vector::one*radius_, points, faces);
    triSurface tri(triFaceList(faces), points);
    refine3D(tri);
    return autoPtr<standAlonePatch>
    (
        new standAlonePatch(tri.faces(), tri.points())
    );
}


bool Foam::immersedSphere::adjustPoints(pointField& points) const
{
    forAll(points, i)
    {
        scalar& x = points[i][0];
        scalar& y = points[i][1];
        scalar& z = points[i][2];

        scalar phi = acos(z/max(mag(points[i]), small));
        scalar theta = atan2(y, x);

        x = radius_*sin(phi)*cos(theta);
        y = radius_*sin(phi)*sin(theta);
        z = radius_*cos(phi);
    }
    return true;
}


Foam::labelList
Foam::immersedSphere::calcInside(const pointField& points) const
{
    scalarList R(mag(points - object_.centre()));
    labelList insidePoints(points.size(), -1);
    label pi = 0;
    forAll(insidePoints, i)
    {
        if (R[i] <= radius_)
        {
            insidePoints[pi++] =  i;
        }
    }
    insidePoints.resize(pi);
    return insidePoints;
}


bool Foam::immersedSphere::inside(const point& pt) const
{
    return (mag(pt - object_.centre()) <= radius_);
}


void Foam::immersedSphere::write(Ostream& os) const
{
    writeEntry(os, "radius", radius_);
}

// ************************************************************************* //
