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

#include "immersedCircle.H"
#include "addToRunTimeSelectionTable.H"
#include "immersedBoundaryObject.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedCircle, 0);
    addToRunTimeSelectionTable(immersedShape, immersedCircle, twoD);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Make list of points axisymmetric and create face indexing
void Foam::immersedCircle::get2DPoints()
{
    scalar pi = Foam::constant::mathematical::pi;

    scalar fullRotation = pi;
    if (this->ai_ == -1)
    {
        fullRotation *= 2.0;
    }
    label nPoints = (fullRotation*radius_)/dx_ + 1;
    scalar dTheta = fullRotation/scalar(nPoints);

    points0_.resize(nPoints + 1);
    forAll(points0_, i)
    {
        scalar theta = scalar(i)*dTheta;
        scalar x = radius_*cos(theta);
        scalar y = radius_*sin(theta);

        points0_[i][xi_] = x;
        points0_[i][yi_] = y;
        points0_[i][ei_] = 0;
    }
}

// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::immersedCircle::immersedCircle
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

    List<face> faces;
    if (ei_ != -1)
    {
        get2DPoints();
        if (ai_ != -1)
        {
            extrudeAxi(points0_, faces);
        }
        else
        {
            extrude2D(points0_, faces);
        }
    }
    patchPtr_.set(new standAlonePatch(faces, points0_));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedCircle::~immersedCircle()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList
Foam::immersedCircle::calcInside(const pointField& points) const
{
    labelList insidePoints(points.size(), -1);
    scalarList R(mag(points - object_.centreOfRotation()));

    label pi = 0;
    forAll(insidePoints, i)
    {
        if (R[i] <= radius_)
        {
            insidePoints[pi++] = i;
        }
    }
    insidePoints.resize(pi);
    return insidePoints;
}



bool Foam::immersedCircle::inside(const point& pt) const
{
    scalar R(mag(pt - object_.centreOfRotation()));
    if (R <= radius_)
    {
        return true;
    }
    return false;
}


// ************************************************************************* //
