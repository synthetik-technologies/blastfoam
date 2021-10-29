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
        dynamicBlastFvMesh::New
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
    fullTriMeshPtr_(nullptr)
{
    read(dict);

    const polyBoundaryMesh& bMesh = fvMesh_.boundaryMesh();
    List<wordRe> names;
    List<wordRe> allNames(bMesh.names());
    if (dict.found("patchName"))
    {
        names.append(dict.lookup<wordRe>("patchName"));
    }
    else
    {
        names.append(List<wordRe>(bMesh.names()));
    }

    labelHashSet patches(bMesh.patchSet(names));
    labelHashSet allPatches(bMesh.patchSet(allNames));
    forAll(bMesh, patchi)
    {
        const polyPatch& patch = bMesh[patchi];
        if
        (
            isA<emptyPolyPatch>(patch)
         || isA<wedgePolyPatch>(patch)
        )
        {
            patches.unset(patch.index());
        }
        else if (isA<processorPolyPatch>(patch))
        {
            allPatches.unset(patch.index());
            patches.unset(patch.index());
        }
    }

    fullTriMeshPtr_ = new triSurface(triangulate(bMesh, allPatches));
    tssPtr_.set(new triSurfaceSearch(*fullTriMeshPtr_));

    autoPtr<triSurface> tri(new triSurface(triangulate(bMesh, patches)));
    patchPtr_.set
    (
        new standAlonePatch
        (
            move(tri->faces()),
            move(tri->points())
        )
    );

    //- points0_ are not used
    points0_.clear();
    writeVTK();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedFvMesh::~immersedFvMesh()
{
    deleteDemandDrivenData(fullTriMeshPtr_);
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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


void Foam::immersedFvMesh::movePoints()
{
    faceCentresOld_ = patchPtr_->faceCentres();

    List<wordRe> names(1, wordRe(this->object_.patchName()));
    const polyBoundaryMesh& bMesh(fvMesh_.boundaryMesh());
    labelHashSet patches(bMesh.patchSet(names));

    List<wordRe> allNames(bMesh.names());
    labelHashSet allPatches(bMesh.patchSet(allNames));
    forAll(bMesh, patchi)
    {
        const polyPatch& patch = bMesh[patchi];
        if (isA<processorPolyPatch>(patch))
        {
            allPatches.unset(patch.index());
        }
    }

    fullTriMeshPtr_ = triangulate(bMesh, allPatches).ptr();
    tssPtr_.reset(new triSurfaceSearch(*fullTriMeshPtr_));
    autoPtr<triSurface> tri(triangulate(bMesh, patches));
    patchPtr_() =
        PrimitivePatch<faceList, const pointField&>
        (
            tri->faces(),
            tri->points()
        );

    bb_ = boundBox(fullTriMeshPtr_->points());
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
    pi = 0;


    boolList inside(tssPtr_->calcInside(validPoints));

    forAll(inside, i)
    {
        if(inside[i]) // Flipped orientation
        {
            insidePoints[pi++] = insidePoints[i];
        }
    }
    insidePoints.resize(pi);
    return insidePoints;
}


bool Foam::immersedFvMesh::inside(const point& pt) const
{
    if (bb_.contains(pt))
    {
        pointField pts(1, pt);
        return tssPtr_->calcInside(pts)[0];
    }
    return false;
}
// ************************************************************************* //
