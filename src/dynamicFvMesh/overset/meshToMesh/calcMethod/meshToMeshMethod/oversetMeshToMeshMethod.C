/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2014 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "oversetMeshToMeshMethod.H"
#include "tetOverlapVolume.H"
#include "OFstream.H"
#include "Time.H"
#include "treeBoundBox.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(oversetMeshToMeshMethod, 0);
    defineRunTimeSelectionTable(oversetMeshToMeshMethod, meshComponents);
}

Foam::scalar Foam::oversetMeshToMeshMethod::tolerance_ = 1e-6;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::oversetMeshToMeshMethod::maskCells() const
{
    boundBox intersectBb
    (
        max(src_.bounds().min(), tgt_.bounds().min()),
        min(src_.bounds().max(), tgt_.bounds().max())
    );

    intersectBb.inflate(0.01);

    const cellList& srcCells = src_.cells();
    const faceList& srcFaces = src_.faces();
    const pointField& srcPts = src_.points();

    DynamicList<label> cells(src_.nCells());
    forAll(srcCells, srcI)
    {
        boundBox cellBb(srcCells[srcI].points(srcFaces, srcPts), false);
        if (intersectBb.overlaps(cellBb))
        {
            cells.append(srcI);
        }
    }

    if (debug)
    {
        Pout<< "participating source mesh cells: " << cells.size() << endl;
    }

    return move(cells);
}


bool Foam::oversetMeshToMeshMethod::intersect
(
    const label srcCelli,
    const label tgtCelli
) const
{
    scalar threshold = tolerance_*src_.cellVolumes()[srcCelli];

    tetOverlapVolume overlapEngine;

    // Note: avoid demand-driven construction of cellPoints
    // treeBoundBox bbTgtCell(tgt_.points(), tgt_.cellPoints()[tgtCelli]);
    const UList<label>& cellFaces = tgt_.cells()[tgtCelli];
    treeBoundBox bbTgtCell(tgt_.points(), tgt_.faces()[cellFaces[0]]);
    for (label i = 1; i < cellFaces.size(); ++i)
    {
        boundBox tmpBb(tgt_.points(), tgt_.faces()[cellFaces[i]], false);
        bbTgtCell.min() = min(bbTgtCell.min(), tmpBb.min());
        bbTgtCell.max() = max(bbTgtCell.max(), tmpBb.max());
    }

    return overlapEngine.cellCellOverlapMinDecomp
    (
        src_,
        srcCelli,
        tgt_,
        tgtCelli,
        bbTgtCell,
        threshold
    );
}


Foam::scalar Foam::oversetMeshToMeshMethod::interVol
(
    const label srcCelli,
    const label tgtCelli
) const
{
    tetOverlapVolume overlapEngine;

    // Note: avoid demand-driven construction of cellPoints
    // treeBoundBox bbTgtCell(tgt_.points(), tgt_.cellPoints()[tgtCelli]);
    const UList<label>& cellFaces = tgt_.cells()[tgtCelli];
    treeBoundBox bbTgtCell(tgt_.points(), tgt_.faces()[cellFaces[0]]);
    for (label i = 1; i < cellFaces.size(); ++i)
    {
        boundBox tmpBb(tgt_.points(), tgt_.faces()[cellFaces[i]], false);
        bbTgtCell.min() = min(bbTgtCell.min(), tmpBb.min());
        bbTgtCell.max() = min(bbTgtCell.max(), tmpBb.max());
    }

    scalar vol = overlapEngine.cellCellOverlapVolumeMinDecomp
    (
        src_,
        srcCelli,
        tgt_,
        tgtCelli,
        bbTgtCell
    );

    return vol;
}


Foam::Tuple2<Foam::scalar, Foam::point>
Foam::oversetMeshToMeshMethod::interVolAndCentroid
(
    const label srcCelli,
    const label tgtCelli
)
{
    NotImplemented;
    return Tuple2<scalar, point>(0.0, Zero);

//     tetOverlapVolume overlapEngine;
//
//     // Note: avoid demand-driven construction of cellPoints
//     // treeBoundBox bbTgtCell(tgt_.points(), tgt_.cellPoints()[tgtCelli]);
//     const UList<label>& cellFaces = tgt_.cells()[tgtCelli];
//     treeBoundBox bbTgtCell(tgt_.points(), tgt_.faces()[cellFaces[0]]);
//     for (label i = 1; i < cellFaces.size(); ++i)
//     {
//         boundBox tmpBb(tgt_.points(), tgt_.faces()[cellFaces[i]]);
//         bbTgtCell.min() = min(bbTgtCell.min(), tmpBb.min());
//         bbTgtCell.max() = min(bbTgtCell.max(), tmpBb.max());
//     }
//
//     Tuple2<scalar, point> volAndInertia =
//     overlapEngine.cellCellOverlapVolumeMinDecomp
//     (
//         src_,
//         srcCelli,
//         tgt_,
//         tgtCelli,
//         bbTgtCell
//     );
//
//     // Convert from inertia to centroid
//     if (volAndInertia.first() <= ROOTVSMALL)
//     {
//         volAndInertia.first() = 0.0;
//         volAndInertia.second() = Zero;
//     }
//     else
//     {
//         volAndInertia.second() /= volAndInertia.first();
//     }
//
//     return volAndInertia;
}


void Foam::oversetMeshToMeshMethod::appendNbrCells
(
    const label celli,
    const polyMesh& mesh,
    const DynamicList<label>& visitedCells,
    DynamicList<label>& nbrCellIDs
) const
{
    const labelList& nbrCells = mesh.cellCells()[celli];

    // filter out cells already visited from cell neighbours
    forAll(nbrCells, i)
    {
        label nbrCelli = nbrCells[i];

        if
        (
            (findIndex(visitedCells, nbrCelli) == -1)
         && (findIndex(nbrCellIDs, nbrCelli) == -1)
        )
        {
            nbrCellIDs.append(nbrCelli);
        }
    }
}


bool Foam::oversetMeshToMeshMethod::initialise
(
    labelListList& srcToTgtAddr,
    scalarListList& srcToTgtWght,
    labelListList& tgtToSrcAddr,
    scalarListList& tgtToSrcWght
) const
{
    srcToTgtAddr.setSize(src_.nCells());
    srcToTgtWght.setSize(src_.nCells());
    tgtToSrcAddr.setSize(tgt_.nCells());
    tgtToSrcWght.setSize(tgt_.nCells());

    if (!src_.nCells())
    {
        return false;
    }
    else if (!tgt_.nCells())
    {
        if (debug)
        {
            Pout<< "mesh interpolation: have " << src_.nCells() << " source "
                << " cells but no target cells" << endl;
        }

        return false;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oversetMeshToMeshMethod::oversetMeshToMeshMethod
(
    const polyMesh& src,
    const polyMesh& tgt
)
:
    src_(src),
    tgt_(tgt),
    V_(0.0)
{
    if (!src_.nCells() || !tgt_.nCells())
    {
        if (debug)
        {
            Pout<< "mesh interpolation: cells not on processor: Source cells = "
                << src_.nCells() << ", target cells = " << tgt_.nCells()
                << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::oversetMeshToMeshMethod::~oversetMeshToMeshMethod()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::oversetMeshToMeshMethod::writeConnectivity
(
    const polyMesh& mesh1,
    const polyMesh& mesh2,
    const labelListList& mesh1ToMesh2Addr
) const
{
    Pout<< "Source size = " << mesh1.nCells() << endl;
    Pout<< "Target size = " << mesh2.nCells() << endl;

    word fName("addressing_" + mesh1.name() + "_to_" + mesh2.name());

    if (Pstream::parRun())
    {
        fName = fName +  "_proc" + Foam::name(Pstream::myProcNo());
    }

    OFstream os(src_.time().path()/fName + ".obj");

    label vertI = 0;
    forAll(mesh1ToMesh2Addr, i)
    {
        const labelList& addr = mesh1ToMesh2Addr[i];
        forAll(addr, j)
        {
            label celli = addr[j];
            const vector& c0 = mesh1.cellCentres()[i];

            const cell& c = mesh2.cells()[celli];
            const pointField pts(c.points(mesh2.faces(), mesh2.points()));
            forAll(pts, j)
            {
                const point& p = pts[j];
                os  << "v " << p.x() << ' ' << p.y() << ' ' << p.z() << nl;
                vertI++;
                os  << "v " << c0.x() << ' ' << c0.y() << ' ' << c0.z()
                    << nl;
                vertI++;
                os  << "l " << vertI - 1 << ' ' << vertI << nl;
            }
        }
    }
}


// ************************************************************************* //
