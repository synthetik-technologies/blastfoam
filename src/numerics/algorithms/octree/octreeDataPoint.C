/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "octreeDataPoint.H"

#include "labelList.H"
#include "treeBoundBox.H"
#include "octree.H"
#include "linePointRef.H"
#include "pointHit.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::octreeDataPoint::octreeDataPoint(const pointField& points)
:
    points_(points)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Get type of volume
Foam::label Foam::octreeDataPoint::getSampleType
(
    const octree<octreeDataPoint>&,
    const point&
) const
{
    return octree<octreeDataPoint>::UNKNOWN;
}


bool Foam::octreeDataPoint::overlaps
(
    const label index,
    const treeBoundBox& sampleBb
) const
{
    return sampleBb.contains(points_[index]);
}


bool Foam::octreeDataPoint::contains
(
    const label,
    const point&
) const
{
    notImplemented
    (
        "octreeDataPoint::contains(const label, const point&)"
    );

    return false;
}


bool Foam::octreeDataPoint::intersects
(
    const label,
    const point&,
    const point&,
    point&
) const
{
    notImplemented
    (
        "octreeDataPoint::intersects(const label, const point&,"
        "const point&, point&)"
    );

    return false;
}


bool Foam::octreeDataPoint::findTightest
(
    const label,
    const point&,
    treeBoundBox&
) const
{
    notImplemented
    (
        "octreeDataPoint::findTightest(const label, const point&,"
        "treeBoundBox&)"
    );

    return false;
}


Foam::scalar Foam::octreeDataPoint::calcSign
(
    const label,
    const point&,
    vector& n
) const
{
    n = vector::zero;

    return 1;
}


// Calculate nearest point on/in shapei
inline Foam::scalar Foam::octreeDataPoint::calcNearest
(
    const label index,
    const point& sample,
    point& nearest
) const
{
    nearest = points_[index];
    return magSqr(points_[index] - sample);
}


void Foam::octreeDataPoint::write
(
    Ostream& os,
    const label index
) const
{
    if ((index < 0) || (index > points().size()))
    {
        FatalErrorIn("octreeDataPoint::write(Ostream&, const label)")
            << "Index " << index << " outside 0.." << points().size()
            << abort(FatalError);
    }
    os << ' ' << points()[index];
}


// Calculate nearest point on/in shapei
Foam::scalar Foam::octreeDataPoint::calcNearest
(
    const label index,
    const linePointRef& ln,
    point& linePt,
    point& shapePt
) const
{
    // Nearest point on shape
    shapePt = points_[index];

    // Nearest point on line
    pointHit pHit = ln.nearestDist(shapePt);

    linePt = pHit.rawPoint();

    return pHit.distance();
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::octreeDataPoint& ocPts
)
{
    return os << ocPts.points();
}


// ************************************************************************* //
