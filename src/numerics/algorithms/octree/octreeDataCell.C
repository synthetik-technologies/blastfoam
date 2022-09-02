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

Description

\*---------------------------------------------------------------------------*/

#include "octreeDataCell.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "treeNode.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::octreeDataCell::octreeDataCell
(
    const polyMesh& mesh,
    const labelList& cellLabels,
    const treeBoundBoxList& bbs
)
:
    mesh_(mesh),
    cellLabels_(cellLabels),
    bbs_(bbs)
{}


// Construct from mesh (assumes all cells)
Foam::octreeDataCell::octreeDataCell
(
    const polyMesh& mesh
)
:
    mesh_(mesh),
    cellLabels_(mesh_.nCells()),
    bbs_
    (
        mesh_.nCells(),
        treeBoundBox::invertedBox
    )
{
    // Set one-one indexing
    for (label i=0; i < mesh_.nCells(); i++)
    {
        cellLabels_[i] = i;
    }

    const pointField& points = mesh_.points();
    const faceList& faces = mesh_.faces();
    const cellList& cells = mesh_.cells();

    forAll(cells, celli)
    {
        const labelList& facesi = cells[celli];

        forAll(facesi, facei)
        {
            const labelList& pointsi = faces[facesi[facei]];

            forAll(pointsi, pointi)
            {
                const point& p = points[pointsi[pointi]];

                bbs_[celli].min() = min(bbs_[celli].min(), p);
                bbs_[celli].max() = max(bbs_[celli].max(), p);
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::octreeDataCell::getSampleType
(
    octree<octreeDataCell>&,
    const point&
) const
{
    return octree<octreeDataCell>::UNKNOWN;
}


bool Foam::octreeDataCell::overlaps
(
    const label index,
    const treeBoundBox& cubeBb
) const
{
    return cubeBb.overlaps(bbs_[index]);
}


bool Foam::octreeDataCell::contains
(
    const label index,
    const point& sample
) const
{
    return mesh_.pointInCell(sample, cellLabels_[index]);
}


bool Foam::octreeDataCell::intersects
(
    const label,
    const point&,
    const point&,
    point&
) const
{
    //Hack: don't know what to do here.

    notImplemented
    (
        "octreeDataCell::intersects(const label, const point&,"
        "const point&, point&)"
    );

    return false;
}


bool Foam::octreeDataCell::findTightest
(
    const label index,
    const point& sample,
    treeBoundBox& tightest
) const
{

    // get nearest and furthest away vertex
    point myNear, myFar;
    bbs_[index].calcExtremities(sample, myNear, myFar);

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
Foam::scalar Foam::octreeDataCell::calcSign
(
    const label,
    const point&,
    vector& n
) const
{
    n = vector::zero;

    return GREAT;
}


// Calculate nearest point on/in shapei
Foam::scalar Foam::octreeDataCell::calcNearest
(
    const label index,
    const point& sample,
    point& nearest
) const
{
    nearest = mesh_.cellCentres()[cellLabels_[index]];

    return mag(nearest - sample);
}


// Calculate nearest point on/in shapei
Foam::scalar Foam::octreeDataCell::calcNearest
(
    const label index,
    const linePointRef& ln,
    point& linePt,
    point& shapePt
) const
{
    notImplemented
    (
        "octreeDataCell::calcNearest(const label, const linePointRef&"
        ", point& linePt, point&)"
    );
    return GREAT;
}


void Foam::octreeDataCell::write
(
    Ostream& os,
    const label index
) const
{
    os << cellLabels_[index] << " " << bbs_[index];
}


// ************************************************************************* //
