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

#include "immersedEllipsoid.H"
#include "addToRunTimeSelectionTable.H"
#include "immersedBoundaryObject.H"
#include "constants.H"
#include "immersedBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedEllipsoid, 0);
    addToRunTimeSelectionTable(immersedShape, immersedEllipsoid, threeD);
}

// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::immersedEllipsoid::immersedEllipsoid
(
    const polyMesh& pMesh,
    const immersedBoundaryObject& ibo,
    const dictionary& dict
)
:
    immersedShape(pMesh, ibo, dict),
    a_(dict.lookup<scalar>("a")),
    b_(dict.lookup<scalar>("b")),
    c_(dict.lookup<scalar>("c"))
{
    read(dict);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedEllipsoid::~immersedEllipsoid()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::standAlonePatch>
Foam::immersedEllipsoid::createPatch() const
{
    List<face> faces;
    pointField points;
    immersedBox::makeBox(vector(a_, b_, c_), points, faces);
    triSurface tri(triSurface(triFaceList(faces), points));

    refine3D(tri);

    return autoPtr<standAlonePatch>
    (
        new standAlonePatch(tri.faces(), tri.points())
    );
}


bool Foam::immersedEllipsoid::adjustPoints(pointField& points) const
{
    forAll(points, i)
    {
        scalar& x = points[i][0];
        scalar& y = points[i][1];
        scalar& z = points[i][2];

        scalar theta = atan2(y, x);
        scalar phi = acos(z/max(mag(points[i]), small));

        scalar r =
            a_*b_*c_
           /sqrt
            (
                sqr(b_*c_*sin(phi)*cos(theta))
              + sqr(a_)
               *(
                    sqr(c_*sin(phi)*sin(theta))
                  + sqr(b_*cos(phi))
                )
            );
        x = r*sin(phi)*cos(theta);
        y = r*sin(phi)*sin(theta);
        z = r*cos(phi);
    }
    return true;
}


Foam::labelList
Foam::immersedEllipsoid::calcInside(const pointField& points) const
{
    labelList insidePoints(points.size(), -1);
    pointField validPoints(points.size());
    label pi = 0;
    forAll(points, i)
    {
        if (bb_.contains(points[i]))
        {
            validPoints[pi] = object_.inverseTransform(points[i]);
            insidePoints[pi++] = i;
        }
    }
    if (!pi)
    {
        return labelList();
    }

    validPoints.resize(pi);
    insidePoints.resize(pi);
    pi = 0;

    forAll(insidePoints, i)
    {
        const scalar& x = validPoints[i][0];
        const scalar& y = validPoints[i][1];
        const scalar& z = validPoints[i][2];
        scalar phi = acos(z/max(mag(validPoints[i]), small));
        scalar theta = atan2(y, x);
        scalar r =
            a_*b_*c_
           /sqrt
            (
                sqr(b_*c_*sin(phi)*cos(theta))
              + sqr(a_)
               *(
                    sqr(c_*sin(phi)*sin(theta))
                  + sqr(b_*cos(phi))
                )
            );
        if (mag(validPoints[i]) <= r)
        {
            insidePoints[pi++] = insidePoints[i];
        }
    }
    if (!pi)
    {
        return labelList();
    }

    insidePoints.resize(pi);
    return insidePoints;
}


bool Foam::immersedEllipsoid::inside(const point& pt) const
{
    if (bb_.contains(pt))
    {

        vector validPt(object_.inverseTransform(pt));

        scalar phi = acos(validPt.z()/max(mag(validPt), small));
        scalar theta = atan2(validPt.y(), validPt.x());
        scalar r =
            a_*b_*c_
           /sqrt
            (
                sqr(b_*c_*sin(phi)*cos(theta))
              + sqr(a_)
               *(
                    sqr(c_*sin(phi)*sin(theta))
                  + sqr(b_*cos(phi))
                )
            );
        if (mag(validPt) <= r)
        {
            return true;
        }
    }
    return false;
}


void Foam::immersedEllipsoid::write(Ostream& os) const
{
    writeEntry(os, "a", a_);
    writeEntry(os, "b", b_);
    writeEntry(os, "c", c_);
}


// ************************************************************************* //
