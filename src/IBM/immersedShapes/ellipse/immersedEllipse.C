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

#include "immersedEllipse.H"
#include "addToRunTimeSelectionTable.H"
#include "immersedBoundaryObject.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedEllipse, 0);
    addToRunTimeSelectionTable(immersedShape, immersedEllipse, twoD);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Make list of points axisymmetric and create face indexing
void Foam::immersedEllipse::get2DPoints()
{
    scalar pi = Foam::constant::mathematical::pi;

    scalar fullRotation = pi;
    if (this->ai_ == -1)
    {
        fullRotation *= 2.0;
    }
    scalar h = sqr(a_ - b_)/sqr(a_ + b_);
    scalar L = (a_ + b_)*(1.0 + 2.0*h/(10.0 + sqrt(4.0 - 3.0*h)))/2.0;

    label nPoints = fullRotation*L/dx_ + 1;
    scalar dTheta = fullRotation/scalar(nPoints);

    points0_.resize(nPoints + 1);
    forAll(points0_, i)
    {
        scalar theta = scalar(i)*dTheta;
        scalar x = a_*cos(theta);
        scalar y = b_*sin(theta);

        points0_[i][xi_] = x;
        points0_[i][yi_] = y;
        points0_[i][ei_] = 0;
    }
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::immersedEllipse::immersedEllipse
(
    const polyMesh& pMesh,
    const immersedBoundaryObject& ibo,
    const dictionary& dict
)
:
    immersedShape(pMesh, ibo, dict),
    a_(dict.lookup<scalar>("a")),
    b_(dict.lookup<scalar>("b"))
{
    read(dict);

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
    patchPtr_.set(new standAlonePatch(faces, points0_));
    correctCentreOfMass();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedEllipse::~immersedEllipse()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::labelList
Foam::immersedEllipse::calcInside(const pointField& points) const
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

    forAll(insidePoints, i)
    {
        scalar theta = atan2(validPoints[i][yi_], validPoints[i][xi_]);
        scalar r =  a_*b_/sqrt(sqr(b_*cos(theta)) + sqr(a_*sin(theta)));

        if (mag(zeroDir(validPoints[i])) <= r)
        {
            insidePoints[pi++] = insidePoints[i];
        }
    }
    insidePoints.resize(pi);
    return insidePoints;
}


bool Foam::immersedEllipse::inside(const point& pt) const
{
    if (bb_.contains(pt))
    {
        point validPt(object_.inverseTransform(pt));
        scalar theta = atan2(validPt[yi_], validPt[xi_]);
        scalar r =  a_*b_/sqrt(sqr(b_*cos(theta)) + sqr(a_*sin(theta)));

        if (mag(zeroDir(validPt)) <= r)
        {
            return true;
        }
    }
    return false;
}


// ************************************************************************* //
