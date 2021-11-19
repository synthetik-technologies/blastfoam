/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "standAlonePatch.H"
#include "vtkWritePolyData.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::standAlonePatch::standAlonePatch
(
    const faceList& faces,
    const pointField& points
)
:
    PrimitivePatch<faceList, pointField>(faces, points),
    pointsRef_(const_cast<pointField&>(this->points()))
{}


Foam::standAlonePatch::standAlonePatch
(
    faceList&& faces,
    pointField&& points
)
:
    PrimitivePatch<faceList, pointField>(faces, points),
    pointsRef_(const_cast<pointField&>(this->points()))
{}


Foam::standAlonePatch::standAlonePatch
(
    faceList&& faces,
    List<point>&& points
)
:
    PrimitivePatch<faceList, pointField>(move(faces), move(points)),
    pointsRef_(const_cast<pointField&>(this->points()))
{}


Foam::standAlonePatch::standAlonePatch
(
    faceList& faces,
    pointField& points,
    const bool reuse
)
:
    PrimitivePatch<faceList, pointField>(faces, points, reuse),
    pointsRef_(const_cast<pointField&>(this->points()))
{}


Foam::standAlonePatch::standAlonePatch(const standAlonePatch& pp)
:
    PrimitivePatch<faceList, pointField>(pp),
    pointsRef_(const_cast<pointField&>(this->points()))
{}


Foam::standAlonePatch::standAlonePatch(standAlonePatch&& pp)
:
    PrimitivePatch<faceList, pointField>(pp),
    pointsRef_(const_cast<pointField&>(this->points()))
{}


Foam::standAlonePatch::standAlonePatch
(
    Istream& is,
    const pointField& points
)
:
    PrimitivePatch<faceList, pointField>(is, points),
    pointsRef_(const_cast<pointField&>(this->points()))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::standAlonePatch::~standAlonePatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::standAlonePatch::movePoints
(
    const pointField& pts
)
{
    if (debug)
    {
        Pout<< "standAlonePatch::"
            << "movePoints() : "
            << "recalculating PrimitivePatch geometry following mesh motion"
            << endl;
    }

    clearGeom();
    pointsRef_ = pts;
}


void Foam::standAlonePatch::writeVTK(const word& name) const
{

    Info<<"writing "<< name << ".vtk" << endl;
    vtkWritePolyData::write
    (
        name + ".vtk",
        name,
        false,
        points(),
        labelList(),
        edgeList(),
        *this
    );
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::standAlonePatch::operator=(const standAlonePatch& pp)
{
    clearOut();

    pointsRef_ = pp.points();
    faceList::operator=(pp);
}


void Foam::standAlonePatch::operator=(standAlonePatch&& pp)
{
    clearOut();

    faceList::operator=(move(pp));
    pointsRef_ = move(pp.pointsRef_);
}

// ************************************************************************* //
