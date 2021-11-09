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

Foam::tmp<Foam::pointField> Foam::immersedCircle::get2DPoints() const
{
    scalar pi = Foam::constant::mathematical::pi;

    scalar fullRotation = pi;
    if (this->ai_ == -1)
    {
        fullRotation *= 2.0;
    }
    label nPoints = (fullRotation*radius_)/dx_ + 1;
    scalar dTheta = fullRotation/scalar(nPoints);

    tmp<pointField> tpoints(new pointField(nPoints + 1));
    pointField& points = tpoints.ref();
    forAll(points, i)
    {
        scalar theta = scalar(i)*dTheta;
        scalar x = radius_*cos(theta);
        scalar y = radius_*sin(theta);

        points[i][xi_] = x;
        points[i][yi_] = y;
        points[i][ei_] = 0;
    }
    return tpoints;
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
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedCircle::~immersedCircle()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::standAlonePatch>
Foam::immersedCircle::createPatch() const
{
    faceList faces;
    pointField points(get2DPoints());
    if (ai_ != -1)
    {
        extrudeAxi(points, faces);
    }
    else
    {
        extrude2D(points, faces);
    }

    return autoPtr<standAlonePatch>
    (
        new standAlonePatch(faces, points)
    );
}


Foam::labelList
Foam::immersedCircle::calcInside(const pointField& points) const
{
    labelList insidePoints(points.size(), -1);

    label pi = 0;
    forAll(insidePoints, i)
    {
        if (mag(zeroDir(points[i] - object_.centre())) <= radius_)
        {
            insidePoints[pi++] = i;
        }
    }
    insidePoints.resize(pi);
    return insidePoints;
}



bool Foam::immersedCircle::inside(const point& pt) const
{
    scalar R(mag(zeroDir(pt - object_.centre())));
    if (R <= radius_)
    {
        return true;
    }
    return false;
}


void Foam::immersedCircle::write(Ostream& os) const
{
    writeEntry(os, "radius", radius_);
}


// ************************************************************************* //
