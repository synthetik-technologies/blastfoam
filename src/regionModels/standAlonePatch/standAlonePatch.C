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


void Foam::standAlonePatch::writeVTK(const fileName& name) const
{
    const fileName file(name + ".vtk");

    Info<<"writing "<< file.name() << endl;
    // Open the file
    std::ofstream os(file, std::ios::binary);

    // Write header
    os  << "# vtk DataFile Version 2.0" << nl
        << name.name(true) << nl
        << "ASCII" << nl
        << "DATASET POLYDATA"
        << nl;

    const pointField& ps = points();

    os  << "POINTS " << ps.size() << " float" << nl;

    // Write vertex coords
    forAll(ps, pointi)
    {
        if (pointi > 0 && (pointi % 10) == 0)
        {
            os  << nl;
        }
        else
        {
            os  << ' ';
        }
        os  << ps[pointi].x() << ' '
            << ps[pointi].y() << ' '
            << ps[pointi].z();
    }
    os  << nl << nl;

    label ne = 0;
    forAll(*this, facei)
    {
        ne += (*this)[facei].size() + 1;
    }

    os  << "POLYGONS " << size() << ' ' << ne << nl;
    forAll(*this, facei)
    {
        os  << (*this)[facei].size() << ' ';
        forAll((*this)[facei], pi)
        {
            os  << operator[](facei)[pi] << ' ';
        }
        os  << nl;
    }
    os  << nl;
}


Foam::standAlonePatch Foam::standAlonePatch::createGlobalPatch() const
{
    if (!Pstream::parRun())
    {
        return *this;
    }

    List<List<point>> gPoints(Pstream::nProcs());
    faceListList gFaces(Pstream::nProcs());

    // Insert my points
    gPoints[Pstream::myProcNo()] = this->points();

    // Insert my faces
    gFaces[Pstream::myProcNo()] = (*this);

    // Communicate points
    Pstream::gatherList(gPoints);
    Pstream::scatterList(gPoints);

    // Communicate faces
    Pstream::gatherList(gFaces);
    Pstream::scatterList(gFaces);

    label nPoints = 0;
    label nFaces = 0;
    forAll(gPoints, proci)
    {
        nPoints += gPoints[proci].size();
        nFaces += gFaces[proci].size();
    }
    pointField ps(nPoints);
    faceList fs(nFaces);

    label pi = 0;
    label fi = 0;
    forAll(gPoints, proci)
    {
        const label start = pi;
        forAll(gPoints[proci], pj)
        {
            ps[pi++] = gPoints[proci][pj];
        }
        forAll(gFaces[proci], fj)
        {
            fs[fi] = gFaces[proci][fj];
            face& f = fs[fi];
            forAll(f, fpi)
            {
                f[fpi] = gFaces[proci][fj][fpi] + start;
            }
            fi++;
        }
    }
    return standAlonePatch(move(fs), move(ps));
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
