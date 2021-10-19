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

#include "octreeDataBoundBox.H"
#include "labelList.H"
#include "polyMesh.H"
#include "octree.H"
#include "polyPatch.H"
#include "triangleFuncs.H"
#include "linePointRef.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::octreeDataBoundBox, 0);

// const Foam::debug::tolerancesSwitch
// Foam::octreeDataBoundBox::tol
// (
//     "octreeDataBoundBoxTol",
//     1e-6
// );

const Foam::scalar
Foam::octreeDataBoundBox::tol
(
    1e-6
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::octreeDataBoundBox::octreeDataBoundBox
(
    const treeBoundBoxList& bbL
)
:
    allBb_(bbL)
{}


Foam::octreeDataBoundBox::octreeDataBoundBox(const octreeDataBoundBox& shapes)
:
    allBb_(shapes.allBb())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::octreeDataBoundBox::~octreeDataBoundBox()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::octreeDataBoundBox::getSampleType
(
    const octree<octreeDataBoundBox>& oc,
    const point& sample
) const
{
    return octree<octreeDataBoundBox>::UNKNOWN;
}


bool Foam::octreeDataBoundBox::overlaps
(
    const label index,
    const treeBoundBox& sampleBb
) const
{
    return sampleBb.overlaps(allBb_[index]);
}


bool Foam::octreeDataBoundBox::contains(const label, const point&) const
{
    notImplemented
    (
        "octreeDataBoundBox::contains(const label, const point&)"
    );

    return false;
}


bool Foam::octreeDataBoundBox::intersects
(
    const label index,
    const point& start,
    const point& end,
    point& intersectionPoint
) const
{
    notImplemented
    (
        "octreeDataBoundBox::intersects(const label, const point&)"
    );
    return false;
}


bool Foam::octreeDataBoundBox::findTightest
(
    const label index,
    const point& sample,
    treeBoundBox& tightest
) const
{
    // Get furthest away vertex
    point myNear, myFar;
    allBb_[index].calcExtremities(sample, myNear, myFar);

    const point dist = myFar - sample;
    scalar myFarDist = mag(dist);

    point tightestNear, tightestFar;
    tightest.calcExtremities(sample, tightestNear, tightestFar);

    scalar tightestFarDist = mag(tightestFar - sample);

    if (tightestFarDist < myFarDist)
    {
        // Keep current tightest.
        return false;
    }
    else
    {
        // Construct bb around sample and myFar
        const point dist2(fabs(dist.x()), fabs(dist.y()), fabs(dist.z()));

        tightest.min() = sample - dist2;
        tightest.max() = sample + dist2;

        return true;
    }
}


// Determine numerical value of sign of sample compared to shape at index
Foam::scalar Foam::octreeDataBoundBox::calcSign
(
    const label index,
    const point& sample,
    point& n
) const
{
    n = vector::zero;

    return 1;
}


// Calculate nearest point on/in shapei
Foam::scalar Foam::octreeDataBoundBox::calcNearest
(
    const label index,
    const point& sample,
    point& nearest
) const
{
    notImplemented
    (
        "octreeDataBoundBox::calcNearest"
        "(const label index, const point& sample, point& nearest)"
    );

    return GREAT;
}


// Calculate nearest point on/in shapei
Foam::scalar Foam::octreeDataBoundBox::calcNearest
(
    const label index,
    const linePointRef& ln,
    point& linePt,
    point& shapePt
) const
{
    notImplemented
    (
        "octreeDataBoundBox::calcNearest"
        "(const label, const linePointRef&, point&, point&)"
    );
    return GREAT;
}


void Foam::octreeDataBoundBox::write(Ostream& os, const label index) const
{
    os << allBb_[index];
}


// ************************************************************************* //
