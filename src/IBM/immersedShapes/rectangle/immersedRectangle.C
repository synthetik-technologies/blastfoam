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

#include "immersedRectangle.H"
#include "addToRunTimeSelectionTable.H"
#include "immersedBoundaryObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedRectangle, 0);
    addToRunTimeSelectionTable(immersedShape, immersedRectangle, twoD);
}

// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::immersedRectangle::immersedRectangle
(
    const polyMesh& pMesh,
    const immersedBoundaryObject& ibo,
    const dictionary& dict
)
:
    immersedShape(pMesh, ibo, dict),
    L_(dict.lookup<scalar>("L")),
    H_(dict.lookup<scalar>("H")),
    origBb_(vector::zero, vector::one)
{
    read(dict);

    vector p1(Zero);
    vector p2(Zero);
    vector p3(Zero);
    vector p4(Zero);

    p1[xi_] = -L_/2.0;
    p1[yi_] = -H_/2.0;

    p2[xi_] = L_/2.0;
    p2[yi_] = -H_/2.0;

    p3[xi_] = L_/2.0;
    p3[yi_] = H_/2.0;

    p4[xi_] = -L_/2.0;
    p4[yi_] = H_/2.0;
    points_ = List<point>({p1, p2, p3, p4});

    boundBox bb(pMesh_.points());
    scalar minE(returnReduce(bb.min()[ei_], minOp<scalar>()));
    scalar maxE(returnReduce(bb.max()[ei_], maxOp<scalar>()));

    origBb_.max() = p3;
    origBb_.max()[ei_] = maxE;
    origBb_.min() = p1;
    origBb_.min()[ei_] = minE;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedRectangle::~immersedRectangle()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::standAlonePatch>
Foam::immersedRectangle::createPatch() const
{
    pointField points;
    points.append(discretizeLine(points_[0], points_[1]));
    points.append(discretizeLine(points_[1], points_[2]));
    points.append(discretizeLine(points_[2], points_[3]));
    points.append(discretizeLine(points_[3], points_[0]));

    List<face> faces;
    if (ai_ != -1)
    {
        extrudeAxi(points, faces);
    }
    else
    {
        extrude2D(points, faces);
    }
    return autoPtr<standAlonePatch>(new standAlonePatch(faces, points));
}


Foam::labelList
Foam::immersedRectangle::calcInside(const pointField& points) const
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
    validPoints = object_.inverseTransform(validPoints);

    forAll(validPoints, i)
    {
        if (origBb_.contains(validPoints[i]))
        {
            insidePoints[pi++] = insidePoints[i];
        }
    }
    insidePoints.resize(pi);
    return insidePoints;
}


bool Foam::immersedRectangle::inside(const point& pt) const
{
    if (bb_.contains(pt))
    {
        if (origBb_.contains(object_.inverseTransform(pt)))
        {
            return true;
        }
    }
    return false;
}


void Foam::immersedRectangle::write(Ostream& os) const
{
    writeEntry(os, "L", L_);
    writeEntry(os, "H", H_);
}


// ************************************************************************* //
