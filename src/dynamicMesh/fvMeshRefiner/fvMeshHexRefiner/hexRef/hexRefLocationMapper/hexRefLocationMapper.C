/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenFOAM Foundation
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

#include "hexRefLocationMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hexRefLocationMapper::hexRefLocationMapper(const polyMesh& mesh)
:
    mesh_(mesh),
    constructMap_(false),

    edgeSplits_(0),
    newEdgeIndices_(0),

    faceSplits_(0),
    newFaceIndices_(0),

    cellSplits_(0),
    newCellIndices_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::hexRefLocationMapper::~hexRefLocationMapper()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hexRefLocationMapper::addSplitEdges
(
    const labelList& edgeMidPoints,
    const labelList& newEdgePoints
)
{
    if (!constructMap_)
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

    labelList splitEdges(edgeMidPoints.size(), -1);
    labelList newPoints(newEdgePoints.size(), -1);
    label ei = 0;
    forAll(splitEdges, edgei)
    {
        if (edgeMidPoints[edgei] >= 0)
        {
            splitEdges[ei] = edgei;
            newPoints[ei++] = newEdgePoints[edgei];
        }
    }
    splitEdges.resize(ei);
    newPoints.resize(ei);

    edgeSplits_.append(splitEdges);
    newEdgeIndices_.append(newPoints);
}


void Foam::hexRefLocationMapper::addSplitFaces
(
    const labelList& faceMidPoints,
    const labelList& newFacePoints
)
{
    if (!constructMap_)
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

    labelList splitFaces(faceMidPoints.size(), -1);
    labelList newPoints(newFacePoints.size(), -1);
    label fi = 0;
    forAll(splitFaces, facei)
    {
        if (faceMidPoints[facei] >= 0)
        {
            splitFaces[fi] = facei;
            newPoints[fi++] = newFacePoints[facei];
        }
    }
    splitFaces.resize(fi);
    newPoints.resize(fi);

    faceSplits_.append(splitFaces);
    newFaceIndices_.append(newPoints);
}


void Foam::hexRefLocationMapper::addSplitCells
(
    const labelList& cellPoints,
    const labelList& newCellPoints
)
{
    if (!constructMap_)
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

    labelList splitCells(cellPoints.size(), -1);
    labelList newPoints(newCellPoints.size(), -1);
    label ci = 0;
    forAll(splitCells, celli)
    {
        if (cellPoints[celli] >= 0)
        {
            splitCells[ci] = celli;
            newPoints[ci++] = newCellPoints[celli];
        }
    }
    splitCells.resize(ci);
    newPoints.resize(ci);

    cellSplits_.append(splitCells);
    newCellIndices_.append(newPoints);
}


const Foam::scalarListList&
Foam::hexRefLocationMapper::edgeWeights() const
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
Foam::hexRefLocationMapper::faceWeights() const
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
Foam::hexRefLocationMapper::cellWeights() const
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


void Foam::hexRefLocationMapper::clearOut()
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


void Foam::hexRefLocationMapper::interpolateMidPoints
(
    pointField& points
) const
{
    if (!constructMap_)
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


// ************************************************************************* //
