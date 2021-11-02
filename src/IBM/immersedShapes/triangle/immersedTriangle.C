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

#include "immersedTriangle.H"
#include "addToRunTimeSelectionTable.H"
#include "immersedBoundaryObject.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedTriangle, 0);
    addToRunTimeSelectionTable(immersedShape, immersedTriangle, twoD);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::immersedTriangle::get2DPoints()
{
    points0_.clear();
    if (ai_ != -1)
    {
        bool onAxis1 = mag(vertices_[0][ri_]) < small;
        bool onAxis2 = mag(vertices_[1][ri_]) < small;
        bool onAxis3 = mag(vertices_[2][ri_]) < small;

        if (!onAxis1 || !onAxis2)
        {
            bool addFinal = onAxis2 && onAxis3;
            points0_.append
            (
                discretizeLine(vertices_[0], vertices_[1], addFinal)
            );
        }
        if (!onAxis2 || !onAxis3)
        {
            bool addFinal = onAxis3 && onAxis1;
            points0_.append
            (
                discretizeLine(vertices_[1], vertices_[2], addFinal)
            );
        }
        if (!onAxis3 || !onAxis1)
        {
            bool addFinal = onAxis1 && onAxis2;
            points0_.append
            (
                discretizeLine(vertices_[2], vertices_[0], addFinal)
            );
        }
        if (onAxis1)
        {
            vertices_[0][ri_] -= 1e-10;
        }
        if (onAxis2)
        {
            vertices_[1][ri_] -= 1e-10;
        }
        if (onAxis3)
        {
            vertices_[2][ri_] -= 1e-10;
        }
    }
    else
    {
        points0_.append(discretizeLine(vertices_[0], vertices_[1]));
        points0_.append(discretizeLine(vertices_[1], vertices_[2]));
        points0_.append(discretizeLine(vertices_[2], vertices_[0]));
    }
    sortPointsPolar(points0_);
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::immersedTriangle::immersedTriangle
(
    const polyMesh& pMesh,
    const immersedBoundaryObject& ibo,
    const dictionary& dict
)
:
    immersedShape(pMesh, ibo, dict),
    vertices_
    (
        {
            dict.lookup<vector>("p1"),
            dict.lookup<vector>("p2"),
            dict.lookup<vector>("p3")
        }
    )
{
    read(dict);

    vector com(Zero);
    forAll(vertices_, i)
    {
        vertices_[i] = zeroDir(vertices_[i]);
        com += vertices_[i];
    }
    com /= scalar(3);
    this->centreOfMass_[xi_] = com[xi_];
    this->centreOfMass_[yi_] = com[yi_];

    List<face> faces;
    get2DPoints();
    if (ai_ != -1)
    {
        extrudeAxi(points0_, faces);
    }
    else
    {
        extrude2D(points0_, faces);
    }
    vertices_ -= centreOfMass_;
    points0_ -= centreOfMass_;

    //- Compute original orientation
    vector x(vertices_[0] - vertices_[1]);
    scalar theta = atan2(x[yi_], x[xi_]);
    this->orientation_ = tensor
    (
        cos(theta), -sin(theta), 0.0,
        sin(theta), cos(theta), 0.0,
        0.0, 0.0, 1.0
    );

    vertices_ = this->orientation_.T() & vertices_;
    points0_ = this->orientation_.T() & points0_;

    patchPtr_.set(new standAlonePatch(faces, points0_));
    correctCentreOfMass();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedTriangle::~immersedTriangle()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList
Foam::immersedTriangle::calcInside(const pointField& points) const
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
    validPoints.resize(pi);
    insidePoints.resize(pi);
    if (!pi)
    {
        return insidePoints;
    }
    pi = 0;

    // Rotate and move original vertices
    pointField verts(object_.transform(vertices_));
    point centre(sum(verts)/3.0);
    labelList vjs({1, 2, 0});
    pointField outsidePoints(3, Zero);

    // Calculate a points outside each of the faces
    forAll(verts, vi)
    {
        vector normal(calcNormal(verts[vi], verts[vjs[vi]]));
        point mid(0.5*(verts[vi] + verts[vjs[vi]]));
        outsidePoints[vi] = mid - normal;
        correctNormal
        (
            normal,
            centre,
            outsidePoints[vi],
            verts[vi],
            verts[vjs[vi]]
        );
        outsidePoints[vi] = mid - normal;
    }

    forAll(insidePoints, i)
    {
        label nHits = 0;
        for (label vi = 0; vi < 3; vi++)
        {
            nHits += intersection
                (
                    validPoints[i],
                    outsidePoints[vi],
                    verts[vi],
                    verts[vjs[vi]]
                ) ? 1 : 0;
        }
        if (nHits % 2 == 1)
        {
            insidePoints[pi++] = insidePoints[i];
        }
    }
    insidePoints.resize(pi);

    return insidePoints;
}


bool Foam::immersedTriangle::inside(const point& pt) const
{
    if (bb_.contains(pt))
    {
        // Rotate and move original vertices
        pointField verts(object_.transform(vertices_));
        point centre(sum(verts)/3.0);
        labelList vjs({1, 2, 0});
        pointField outsidePoints(3, Zero);

        // Calculate a points outside each of the faces
        forAll(verts, vi)
        {
            vector normal(calcNormal(verts[vi], verts[vjs[vi]]));
            point mid(0.5*(verts[vi] + verts[vjs[vi]]));
            outsidePoints[vi] = mid - normal;
            correctNormal
            (
                normal,
                centre,
                outsidePoints[vi],
                verts[vi],
                verts[vjs[vi]]
            );
            outsidePoints[vi] = mid - normal;
        }

        label nHits = 0;
        for (label vi = 0; vi < 3; vi++)
        {
            nHits += intersection
                (
                    pt,
                    outsidePoints[vi],
                    verts[vi],
                    verts[vjs[vi]]
                ) ? 1 : 0;
        }
        if (nHits % 2 == 1)
        {
            return true;
        }
    }
    return false;
}


// ************************************************************************* //
