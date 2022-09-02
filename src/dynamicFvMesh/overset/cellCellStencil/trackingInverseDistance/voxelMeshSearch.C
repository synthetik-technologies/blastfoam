/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "voxelMeshSearch.H"
#include "polyMesh.H"
#include "processorPolyPatch.H"
#include "IOobject.H"
#include "fvMesh.H"
#include "block.H"
#include "emptyPolyPatch.H"
#include "fvMeshTools.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(voxelMeshSearch, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelVector Foam::voxelMeshSearch::offset
(
    const labelVector& nDivs
)
{
    return labelVector(1, nDivs.x(), nDivs.x()*nDivs.y());
}


Foam::label Foam::voxelMeshSearch::index
(
    const labelVector& nDivs,
    const labelVector& voxel
)
{
    return voxel.x()+voxel.y()*nDivs.x()+voxel.z()*nDivs.x()*nDivs.y();
}


Foam::labelVector Foam::voxelMeshSearch::index3
(
    const labelVector& nDivs,
    label voxeli
)
{
    const label nxy = nDivs.x()*nDivs.y();

    labelVector voxel;
    voxel.z() = voxeli/nxy;
    voxeli = voxeli % nxy;
    voxel.y() = voxeli/nDivs.x();
    voxel.x() = voxeli%nDivs.x();

    return voxel;
}


Foam::labelVector Foam::voxelMeshSearch::index3
(
    const boundBox& bb,
    const labelVector& g,
    const point& pt
)
{
    const vector s(cmptDivide(bb.span(), vector(g.x(), g.y(), g.z())));

    labelVector v
    (
        floor((pt.x()-bb.min().x())/s.x()),
        floor((pt.y()-bb.min().y())/s.y()),
        floor((pt.z()-bb.min().z())/s.z())
    );

    return v;
}


Foam::label Foam::voxelMeshSearch::index
(
    const boundBox& bb,
    const labelVector& g,
    const point& pt,
    const bool clip
)
{
    const vector s(cmptDivide(bb.span(), vector(g.x(), g.y(), g.z())));

    labelVector v
    (
        floor((pt.x()-bb.min().x())/s.x()),
        floor((pt.y()-bb.min().y())/s.y()),
        floor((pt.z()-bb.min().z())/s.z())
    );

    if (clip)
    {
        v[0] = max(0, min(g[0]-1, v[0]));
        v[1] = max(0, min(g[1]-1, v[1]));
        v[2] = max(0, min(g[2]-1, v[2]));
    }
    else if
    (
        v[0] < 0
     || v[1] < 0
     || v[2] < 0
     || v[0] >= g[0]
     || v[1] >= g[1]
     || v[2] >= g[2]
    )
    {
        return -1;
    }

    return index(g, v);
}


Foam::point Foam::voxelMeshSearch::centre
(
    const boundBox& bb,
    const labelVector& g,
    const labelVector& voxel
)
{
    const vector s(cmptDivide(bb.span(), vector(g.x(), g.y(), g.z())));

    return bb.min()+0.5*s+point(voxel[0]*s[0], voxel[1]*s[1], voxel[2]*s[2]);
}


void Foam::voxelMeshSearch::writeGrid
(
    OBJstream& os,
    const boundBox& bb,
    const labelVector& g
)
{
    const vector s(cmptDivide(bb.span(), vector(g.x(), g.y(), g.z())));

    for (label i = 1; i < g[0]; i++)
    {
        for (label j = 0; j < g[1]; j++)
        {
            for (label k = 0; k < g[2]; k++)
            {
                point p1(bb.min()+point((i-1)*s[0], j*s[1], k*s[2]));
                point p2(bb.min()+point(i*s[0], j*s[1], k*s[2]));
                os.write(linePointRef(p1, p2));
            }
        }
    }
    for (label i = 0; i < g[0]; i++)
    {
        for (label j = 1; j < g[1]; j++)
        {
            for (label k = 0; k < g[2]; k++)
            {
                point p1(bb.min()+point(i*s[0], (j-1)*s[1], k*s[2]));
                point p2(bb.min()+point(i*s[0], j*s[1], k*s[2]));
                os.write(linePointRef(p1, p2));
            }
        }
    }
    for (label i = 0; i < g[0]; i++)
    {
        for (label j = 0; j < g[1]; j++)
        {
            for (label k = 1; k < g[2]; k++)
            {
                point p1(bb.min()+point(i*s[0], j*s[1], (k-1)*s[2]));
                point p2(bb.min()+point(i*s[0], j*s[1], k*s[2]));
                os.write(linePointRef(p1, p2));
            }
        }
    }
}


Foam::label Foam::voxelMeshSearch::searchProcPatch
(
    const label faceID,
    const point& searchPoint
) const
{
    const pointField& cellCentres = mesh_.cellCentres();
    const polyBoundaryMesh& bMeshes = mesh_.boundaryMesh();

    label patchi = bMeshes.patchID()[faceID-mesh_.nInternalFaces()];
    const polyPatch& bMeshPatch = bMeshes[patchi];

    if (!isA<processorPolyPatch>(bMeshPatch))
    {
        return -1;
    }
    else
    {
        // Find nearest cell. Linear search since cheaper than constructing
        // tree?
        const labelUList& faceCells = bMeshPatch.faceCells();
        scalar minProximity = GREAT;

        label nearestCellI = -1;
        forAll(faceCells, i)
        {
            const point& cc = cellCentres[faceCells[i]];
            scalar proximity = magSqr(cc-searchPoint);
            if (proximity < minProximity)
            {
                minProximity = proximity;
                nearestCellI = faceCells[i];
            }
        }
        return nearestCellI;
    }
}


Foam::label Foam::voxelMeshSearch::findIntersectedFace
(
    const label cellI,
    const point& p
) const
{
    // Return -1 or the label of the face intersected when tracking from
    // p to centre of cellI

    const faceList& faces = mesh_.faces();
    const pointField& faceCentres = mesh_.faceCentres();
    const pointField& points = mesh_.points();

    const point& cc = mesh_.cellCentres()[cellI];
    const labelList& cFaces = mesh_.cells()[cellI];

    const vector q(cc-p);

    forAll(cFaces, cFacei)
    {
        label facei = cFaces[cFacei];

        pointHit hitInfo = faces[facei].intersection
        (
            p,
            q,
            faceCentres[facei],
            points,
            intersection::algorithm::halfRay
        );

        if (hitInfo.hit() && (hitInfo.distance() < 1))
        {
            return facei;
        }
    }
    return -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::voxelMeshSearch::voxelMeshSearch
(
    const polyMesh& mesh,
    const bool doUpdate
)
:
    mesh_(mesh)
{
    // Determine number of voxels from number of cells in mesh
    const labelVector& dim = mesh_.geometricD();

    // Guarantee at least one voxel
    label nCells = max(1, mesh_.nCells());

    label nDivs = -1;
    if (mesh_.nGeometricD() == 1)
    {
        nDivs = nCells;
    }
    else if (mesh_.nGeometricD() == 2)
    {
        nDivs = label(Foam::sqrt(scalar(nCells)));
    }
    else
    {
        nDivs = label(Foam::cbrt(scalar(nCells)));
    }

    nDivs_ = labelVector(nDivs, nDivs, nDivs);
    forAll(dim, i)
    {
        if (dim[i] == -1)
        {
            nDivs_[i] = 1;
        }
    }

    // Redo the local bounding box
    localBb_ = boundBox(mesh_.points(), false);

    const point eps(1e-10, 1e-10, 1e-10);

    localBb_.min() = localBb_.min()-eps;
    localBb_.max() = localBb_.max()+eps;

    if (debug)
    {
        Pout<< "voxelMeshSearch : mesh:" << mesh_.name()
            << " nDivs:" << nDivs_ << endl;
    }

    if (doUpdate)
    {
        update();
    }
}


Foam::voxelMeshSearch::voxelMeshSearch
(
    const polyMesh& mesh,
    const boundBox& localBb,
    const labelVector& nDivs,
    const bool doUpdate
)
:
    mesh_(mesh),
    localBb_(localBb),
    nDivs_(nDivs)
{
    if (doUpdate)
    {
        update();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::voxelMeshSearch::update()
{
    // Initialise seed cell array

    seedCell_.setSize(cmptProduct(nDivs_));
    seedCell_ = -1;


    // Find seed cells
    const pointField& points = mesh_.points();
    const labelListList& cellPoints = mesh_.cellPoints();

    forAll(cellPoints, celli)
    {
        const labelList& cPoints = cellPoints[celli];

        // Get cell bounding box
        boundBox bb(points, cPoints, false);

        fill(seedCell_, localBb_, nDivs_, bb, celli);
    }


    if (debug)
    {
        Pout<< "voxelMeshSearch : mesh:" << mesh_.name()
            << " nDivs:" << nDivs_
            << " localBb:" << localBb_ << endl;
    }


    //// Small optimisation: make sure the cell centre at least always
    //// returns the cell itself
    //const pointField& cellCentres = mesh_.cellCentres();
    //forAll(cellCentres, celli)
    //{
    //    label voxeli = index(cellCentres[celli]);
    //    seedCell_[voxeli] = celli;
    //}

    return true;
}


Foam::label Foam::voxelMeshSearch::findCell(const point& p) const
{
    // First check if the point is contained in the bounding box, else exit
    if (!localBb_.contains(p))
    {
        return -1;
    }


    // Locate the voxel index for this point. Do not clip.
    label voxeli = index(localBb_, nDivs_, p, false);

    // The point may still be inside the bb but outside the actual domain.
    if (voxeli < 0)
    {
        return -1;
    }
    else
    {
        // Inverse map to compute the seed cell.
        label celli = seedCell_[voxeli];

        if (celli < 0)
        {
            return -1;
        }
        else
        {
            // Simplified, non-parallel tracking from cell centre of
            // celli to wanted location p. Note that the cell thus
            // found does not have to be the absolute 'correct' one as
            // long as at least one of the processors finds a cell.

            track_.clear();
            while (true)
            {
                if (track_.size() < 5)
                {
                    track_.append(celli);
                }

                // I am in celli now. How many faces do I have ?
                label facei = findIntersectedFace(celli, p);

                if (facei == -1)
                {
                    return celli;
                }

                const label startOfTrack(max(0, track_.size()-5));

                label nextCell;
                if (mesh_.isInternalFace(facei))
                {
                    label own = mesh_.faceOwner()[facei];
                    label nei = mesh_.faceNeighbour()[facei];
                    nextCell = (own == celli ? nei : own);

                    if (findIndex(track_, nextCell, startOfTrack) >= 0)
                    {
                        return celli;
                    }
                }
                else
                {
                    nextCell = searchProcPatch(facei, p);

                    if (nextCell == -1 || nextCell == celli)
                    {
                        return nextCell;
                    }
                    else if (findIndex(track_, nextCell, startOfTrack) >= 0)
                    {
                        return -1;  // point is really out
                    }
                }

                celli = nextCell;
            }
            return -1;
        }
    }
}


Foam::autoPtr<Foam::fvMesh> Foam::voxelMeshSearch::makeMesh
(
    const IOobject& io
) const
{
    const cellModel& hex = *(cellModeller::lookup("hex"));

    PtrList<cellShape> cellShapes;
    faceListList boundary;
    pointField points;
    {
        //Info<< "Creating block" << endl;

        block b
        (
            blockDescriptor
            (
                cellShape(hex, identity(8), false),
                localBb_.points(),
                blockEdgeList(),
                blockFaceList(),
                nDivs_,
                List<gradingDescriptors>(12)
            )
        );

        List<FixedList<label, 8>> bCells(b.cells());
        cellShapes.setSize(bCells.size());
        forAll(cellShapes, celli)
        {
            cellShapes.set
            (
                celli,
                new cellShape(hex, labelList(bCells[celli]), false)
            );
        }

        //Info<< "Creating boundary faces" << endl;

        boundary.setSize(b.boundaryPatches().size());
        forAll(boundary, patchi)
        {
            faceList faces(b.boundaryPatches()[patchi].size());
            forAll(faces, facei)
            {
                faces[facei] = face(labelList(b.boundaryPatches()[patchi][facei]));
            }
            boundary[patchi].transfer(faces);
        }

        points.transfer(const_cast<pointField&>(b.points()));
    }

    //Info<< "Creating patch dictionaries" << endl;
    wordList patchNames(boundary.size());
    forAll(patchNames, patchi)
    {
        patchNames[patchi] =
            patchi < 0
          ? word("patch", false)
          : word("patch" + std::to_string(patchi), false);

    }

    PtrList<dictionary> boundaryDicts(boundary.size());
    forAll(boundaryDicts, patchi)
    {
        boundaryDicts.set(patchi, new dictionary());
        dictionary& patchDict = boundaryDicts[patchi];
        patchDict.add("type", emptyPolyPatch::typeName);
    }

    //Info<< "Creating polyMesh" << endl;
    IOobject polyIO(io);
    polyIO.readOpt() = IOobject::NO_READ;
    polyMesh mesh
    (
        //IOobject
        //(
        //    polyMesh::defaultRegion,
        //    runTime.constant(),
        //    runTime,
        //    IOobject::NO_READ
        //),
        polyIO,
        std::move(points),
        cellShapeList(move(cellShapes)),
        boundary,
        patchNames,
        boundaryDicts,
        "defaultFaces",
        emptyPolyPatch::typeName,
        false
    );

    //Info<< "Writing polyMesh" << endl;
    mesh.write();

    //Info<< "Reading fvMesh" << endl;

    createDummyFvMeshFiles
    (
        io.db(),
        io.name(),
        false
    );
    IOobject fvIO(io);
    fvIO.readOpt() = IOobject::MUST_READ;

    return autoPtr<fvMesh>(new fvMesh(fvIO));
}


void Foam::voxelMeshSearch::createDummyFvMeshFiles
(
    const objectRegistry& mesh,
    const word& regionName,
    const bool verbose
) const
{
    // Create dummy system/fv*
    {
        IOobject io
        (
            "fvSchemes",
            mesh.time().system(),
            regionName,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        );

        if (!io.typeHeaderOk<IOdictionary>(false))
        {
            if (verbose)
            {
                Info<< "Writing dummy " << regionName/io.name() << endl;
            }
            dictionary dummyDict;
            dictionary divDict;
            dummyDict.add("divSchemes", divDict);
            dictionary gradDict;
            dummyDict.add("gradSchemes", gradDict);
            dictionary laplDict;
            dummyDict.add("laplacianSchemes", laplDict);

            IOdictionary(io, dummyDict).regIOobject::write();
        }
    }
    {
        IOobject io
        (
            "fvSolution",
            mesh.time().system(),
            regionName,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        );

        if (!io.typeHeaderOk<IOdictionary>(false))
        {
            if (verbose)
            {
                Info<< "Writing dummy " << regionName/io.name() << endl;
            }
            dictionary dummyDict;
            IOdictionary(io, dummyDict).regIOobject::write();
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
