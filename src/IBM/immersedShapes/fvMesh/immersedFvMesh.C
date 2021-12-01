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

#include "immersedFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "immersedBoundaryObject.H"
#include "wedgePolyPatch.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "globalPolyPatch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedFvMesh, 0);
    addToRunTimeSelectionTable(immersedShape, immersedFvMesh, twoD);
    addToRunTimeSelectionTable(immersedShape, immersedFvMesh, threeD);
}



// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::immersedFvMesh::immersedFvMesh
(
    const polyMesh& pMesh,
    const immersedBoundaryObject& ibo,
    const dictionary& dict
)
:
    immersedShape(pMesh, ibo, dict),
    fvMeshPtr_
    (
        dynamicFvMesh::New
        (
            IOobject
            (
                this->name(),
                pMesh.time().timeName(),
                pMesh.time(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            )
        )
    ),
    fvMesh_(fvMeshPtr_()),
    patchName_(dict.lookupOrDefault<word>("patchName", this->name_))
{
    read(dict);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedFvMesh::~immersedFvMesh()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::standAlonePatch>
Foam::immersedFvMesh::createPatch() const
{
    const polyBoundaryMesh& bMesh = fvMesh_.boundaryMesh();
    globalPolyPatch gpp(dictionary(), bMesh[patchName_]);
    return autoPtr<standAlonePatch>
    (
        new standAlonePatch(gpp.globalPatch())
    );
}


Foam::autoPtr<Foam::triSurface> Foam::immersedFvMesh::triangulate
(
    const polyBoundaryMesh& bMesh,
    const labelHashSet& patches
) const
{
    if (Pstream::parRun())
    {
        triSurface tMesh(triSurfaceTools::triangulate(bMesh, patches));

        List<List<labelledTri>> triFaces(Pstream::nProcs());
        triFaces[Pstream::myProcNo()] = tMesh;
        Pstream::gatherList(triFaces);
        Pstream::scatterList(triFaces);

        List<pointField> triPoints(Pstream::nProcs());
        triPoints[Pstream::myProcNo()] = tMesh.points();
        Pstream::gatherList(triPoints);
        Pstream::scatterList(triPoints);

        List<labelledTri> faces(triFaces[0]);
        pointField points(triPoints[0]);
        label startI = points.size();
        for (label i = 1; i < triFaces.size(); i++)
        {
            List<labelledTri>& subFaces(triFaces[i]);
            forAll(subFaces, facei)
            {
                subFaces[facei][0] += startI;
                subFaces[facei][1] += startI;
                subFaces[facei][2] += startI;
            }
            faces.append(triFaces[i]);
            points.append(triPoints[i]);
            startI += triPoints[i].size();
        }
        points *= 1 + small;
        autoPtr<triSurface> tMeshPtr(new triSurface(faces, points));
        tMeshPtr->cleanup(false);

        return tMeshPtr;
    }
    else
    {
        return autoPtr<triSurface>
        (
            new triSurface(triSurfaceTools::triangulate(bMesh, patches))
        );
    }
}


void Foam::immersedFvMesh::movePoints(const pointField& points)
{
    const polyBoundaryMesh& bMesh = fvMesh_.boundaryMesh();
    globalPolyPatch gpp(dictionary(), bMesh[patchName_]);
    const_cast<pointField&>(points) = gpp.globalPatch().points();
    this->centre_ = sum
        (
            fvMesh_.C().primitiveField()
           *fvMesh_.V().field()
        )/sum(fvMesh_.V().field());

    bb_ = boundBox(fvMesh_.points());
    if (ai_ != -1)
    {
        bb_.min()[ai_] = -great;
        bb_.max()[ai_] = great;
    }
    if (ei_ != -1)
    {
        bb_.min()[ei_] = -great;
        bb_.max()[ei_] = great;
    }
}


Foam::labelList
Foam::immersedFvMesh::calcInside(const pointField& points) const
{
    labelList insidePoints(points.size(), -1);
    pointField validPoints(points);
    label pi = 0;
    forAll(points, i)
    {
        if (bb_.contains(points[i]))
        {
            validPoints[pi] = points[i];
            insidePoints[pi++] = i;
        }
    }
    validPoints.resize(pi);
    insidePoints.resize(pi);
    if (!pi)
    {
        return insidePoints;
    }

    boolList inside(validPoints.size());

    forAll(inside, i)
    {
        inside[i] = fvMesh_.findCell(validPoints[i]) >= 0;
    }
    reduce(inside, sumOp<boolList>());

    pi = 0;
    forAll(inside, i)
    {
        if (inside[i])
        {
            insidePoints[pi++] = insidePoints[i];
        }
    }
    insidePoints.resize(pi);
    return insidePoints;
}


Foam::diagTensor Foam::immersedFvMesh::momentOfInertia() const
{
    NotImplemented;
    return diagTensor(great, great, great);
}


bool Foam::immersedFvMesh::inside(const point& pt) const
{
    if (bb_.contains(pt))
    {
        return returnReduce(fvMesh_.findCell(pt) >= 0, orOp<bool>());
    }
    return false;
}


void Foam::immersedFvMesh::write(Ostream& os) const
{}


// ************************************************************************* //
