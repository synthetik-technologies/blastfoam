/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applies Technologies
-------------------------------------------------------------------------------
License
    This file is a derived work of OpenFOAM.

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

#include "locationMapper.H"
#include "pointFields.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(locationMapper, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::locationMapper::locationMapper(const polyMesh& mesh)
:
    regIOobject
    (
        IOobject
        (
            typeName,
            mesh.time().name(),
            mesh
        )
    ),
    mesh_(mesh),
    needMap_(false),
    interpolatedFields_(),

    edgeSplits_(0),
    newEdgeIndices_(0),

    faceSplits_(0),
    newFaceIndices_(0),

    cellSplits_(0),
    newCellIndices_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::locationMapper::~locationMapper()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::locationMapper::addSplitEdges
(
    const labelList& splitEdges,
    const labelList& newEdgePoints
)
{
    if (!needMap())
    {
        return;
    }

    if (!pointsOldPtr_.valid())
    {
        pointsOldPtr_.set(new pointField(mesh_.points()));
    }
    if (!edgesOldPtr_.valid())
    {
        edgesOldPtr_.set(new edgeList(mesh_.edges()));
    }

    labelList newPoints(splitEdges.size(), -1);
    forAll(splitEdges, ei)
    {
        newPoints[ei] = newEdgePoints[splitEdges[ei]];
    }

    edgeSplits_.append(splitEdges);
    newEdgeIndices_.append(newPoints);
}


void Foam::locationMapper::addSplitFaces
(
    const labelList& splitFaces,
    const labelList& newFacePoints
)
{
    if (!needMap())
    {
        return;
    }

    if (!pointsOldPtr_.valid())
    {
        pointsOldPtr_.set(new pointField(mesh_.points()));
    }
    if (!facesOldPtr_.valid())
    {
        facesOldPtr_.set(new faceList(mesh_.faces()));
    }

    labelList newPoints(splitFaces.size(), -1);
    forAll(splitFaces, fi)
    {
        newPoints[fi] = newFacePoints[splitFaces[fi]];
    }

    faceSplits_.append(splitFaces);
    newFaceIndices_.append(newPoints);
}


void Foam::locationMapper::addSplitCells
(
    const labelList& splitCells,
    const labelList& newCellPoints
)
{
    if (!needMap())
    {
        return;
    }

    if (!pointsOldPtr_.valid())
    {
        pointsOldPtr_.set(new pointField(mesh_.points()));
    }
    if (!facesOldPtr_.valid())
    {
        facesOldPtr_.set(new faceList(mesh_.faces()));
    }
    if (!cellsOldPtr_.valid())
    {
        cellsOldPtr_.set(new cellList(mesh_.cells()));
    }

    labelList newPoints(splitCells.size(), -1);
    forAll(splitCells, ci)
    {
        newPoints[ci] = newCellPoints[splitCells[ci]];
    }

    cellSplits_.append(splitCells);
    newCellIndices_.append(newPoints);
}


const Foam::scalarListList&
Foam::locationMapper::edgeWeights() const
{
    if (edgeWeightsPtr_.valid())
    {
        return edgeWeightsPtr_();
    }

    edgeWeightsPtr_.set(new scalarListList(edgeSplits_.size()));
    scalarListList& edgeWs(edgeWeightsPtr_());

    const pointField& pointsOld(pointsOldPtr_());
    const edgeList& edges = edgesOldPtr_();
    forAll(edgeSplits_, i)
    {
        const edge& e = edges[edgeSplits_[i]];
        point pt = e.centre(pointsOld);
        scalarField w(e.size());
        scalar sumWeight = 0.0;
        forAll(e, pti)
        {
            w[pti] = 1.0/mag(pointsOld[e[pti]] - pt);
            sumWeight += w[pti];
        }
        w /= sumWeight;
        edgeWs[i].transfer(w);
    }
    return edgeWs;
}


const Foam::scalarListList&
Foam::locationMapper::faceWeights() const
{
    if (faceWeightsPtr_.valid())
    {
        return faceWeightsPtr_();
    }

    faceWeightsPtr_.set(new scalarListList(faceSplits_.size()));
    scalarListList& faceWs(faceWeightsPtr_());

    const pointField& pointsOld(pointsOldPtr_());
    const faceList& faces = facesOldPtr_();
    forAll(faceSplits_, i)
    {
        const face& f = faces[faceSplits_[i]];
        point pt = f.centre(pointsOld);
        scalarField w(f.size());
        scalar sumWeight = 0.0;
        forAll(f, pti)
        {
            w[pti] = 1.0/mag(pointsOld[f[pti]] - pt);
            sumWeight += w[pti];
        }
        w /= sumWeight;
        faceWs[i].transfer(w);
    }
    return faceWs;
}


const Foam::scalarListList&
Foam::locationMapper::cellWeights() const
{
    if (cellWeightsPtr_.valid())
    {
        return cellWeightsPtr_();
    }

    cellWeightsPtr_.set(new scalarListList(cellSplits_.size()));
    scalarListList& cellWs(cellWeightsPtr_());

    const pointField& pointsOld(pointsOldPtr_());
    const faceList& faces = facesOldPtr_();
    const cellList& cells = cellsOldPtr_();
    forAll(cellSplits_, i)
    {
        const cell& c = cells[cellSplits_[i]];
        point pt = c.centre(pointsOld, faces);
        const labelList cp = c.labels(faces);
        scalarField w(cp.size());
        scalar sumWeight = 0.0;
        forAll(cp, pti)
        {
            w[pti] = 1.0/mag(pointsOld[cp[pti]] - pt);
            sumWeight += w[pti];
        }
        w /= sumWeight;
        cellWs[i].transfer(w);
    }
    return cellWs;
}


void Foam::locationMapper::clearOut() const
{
    pointsOldPtr_.clear();
    edgesOldPtr_.clear();
    facesOldPtr_.clear();
    cellsOldPtr_.clear();

    edgeSplits_.clear();
    newEdgeIndices_.clear();

    faceSplits_.clear();
    newFaceIndices_.clear();

    cellSplits_.clear();
    newCellIndices_.clear();
}


void Foam::locationMapper::interpolateMidPoints
(
    pointField& points
) const
{
    if (!needMap())
    {
        return;
    }

    pointField pointsOld(points);
    points.resize(mesh_.points().size());

    if (edgesOldPtr_.valid())
    {
        const edgeList& edges = edgesOldPtr_();
        forAll(edgeSplits_, i)
        {
            points[newEdgeIndices_[i]] =
                edges[edgeSplits_[i]].centre(pointsOld);
        }
    }

    if (facesOldPtr_.valid())
    {
        const faceList& faces = facesOldPtr_();
        forAll(faceSplits_, i)
        {
            points[newFaceIndices_[i]] =
                faces[faceSplits_[i]].centre(pointsOld);
        }

        if (cellsOldPtr_.valid())
        {
            const cellList& cells = cellsOldPtr_();
            forAll(cellSplits_, i)
            {
                points[newCellIndices_[i]] =
                    cells[cellSplits_[i]].centre(pointsOld, faces);
            }
        }
    }
}


bool Foam::locationMapper::found(const polyMesh& mesh)
{
    return mesh.foundObject<locationMapper>(locationMapper::typeName);
}


const Foam::locationMapper& Foam::locationMapper::New(const polyMesh& mesh)
{
    if (!mesh.foundObject<locationMapper>(locationMapper::typeName))
    {
        locationMapper* lm = new locationMapper(mesh);
        lm->store(lm);
    }
    return mesh.lookupObject<locationMapper>(locationMapper::typeName);
}


Foam::locationMapper& Foam::locationMapper::NewRef(const polyMesh& mesh)
{
    if (!mesh.foundObject<locationMapper>(locationMapper::typeName))
    {
        locationMapper* lm = new locationMapper(mesh);
        lm->store(lm);
    }
    return mesh.lookupObjectRef<locationMapper>(locationMapper::typeName);
}


Foam::locationMapper& Foam::locationMapper::New(polyMesh& mesh)
{
    if (!mesh.foundObject<locationMapper>(locationMapper::typeName))
    {
        locationMapper* lm = new locationMapper(mesh);
        lm->store(lm);
    }
    return mesh.lookupObjectRef<locationMapper>(locationMapper::typeName);
}


// ************************************************************************* //
