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

#include "immersedStl.H"
#include "addToRunTimeSelectionTable.H"
#include "PackedBoolList.H"
#include "polyMesh.H"
#include "triSurface.H"
#include "featureEdgeMesh.H"
#include "extendedFeatureEdgeMesh.H"
#include "triSurfaceSearch.H"
#include "booleanSurface.H"
#include "orientedSurface.H"
#include "immersedBoundaryObject.H"
#include "intersectedSurface.H"
#include "edgeSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedStl, 0);
    addToRunTimeSelectionTable(immersedShape, immersedStl, twoD);
    addToRunTimeSelectionTable(immersedShape, immersedStl, threeD);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::immersedStl::intersectSurfaces
(
    triSurface& surf,
    edgeIntersections& edgeCuts
)
{
    bool hasMoved = false;

    for (label iter = 0; iter < 10; iter++)
    {
        // Determine surface edge intersections. Allow surface to be moved.

        // Number of iterations needed to resolve degenerates
        label nIters = 0;
        {
            triSurfaceSearch querySurf(surf);

            scalarField surfPointTol
            (
                1e-6*edgeIntersections::minEdgeLength(surf)
            );

            // Determine raw intersections
            edgeCuts = edgeIntersections
            (
                surf,
                querySurf,
                surfPointTol
            );

            // Shuffle a bit to resolve degenerate edge-face hits
            {
                pointField points(surf.points());

                nIters =
                    edgeCuts.removeDegenerates
                    (
                        5,              // max iterations
                        surf,
                        querySurf,
                        surfPointTol,
                        points         // work array
                    );

                if (nIters != 0)
                {
                    // Update geometric quantities
                    surf.movePoints(points);
                    hasMoved = true;
                }
            }
        }
    }

    return hasMoved;
}


// Keep on shuffling surface points until no more degenerate intersections.
// Moves both surfaces and updates set of edge cuts.
bool Foam::immersedStl::intersectSurfaces
(
    triSurface& surf1,
    edgeIntersections& edgeCuts1,
    triSurface& surf2,
    edgeIntersections& edgeCuts2
)
{
    bool hasMoved1 = false;
    bool hasMoved2 = false;

    for (label iter = 0; iter < 20; iter++)
    {
        // Determine surface1 edge intersections. Allow surface to be moved.

        // Number of iterations needed to resolve degenerates
        label nIters1 = 0;
        {
            triSurfaceSearch querySurf2(surf2);

            scalarField surf1PointTol
            (
                1e-3*edgeIntersections::minEdgeLength(surf1)
            );

            // Determine raw intersections
            edgeCuts1 = edgeIntersections
            (
                surf1,
                querySurf2,
                surf1PointTol
            );

            // Shuffle a bit to resolve degenerate edge-face hits
            {
                pointField points1(surf1.points());

                nIters1 =
                    edgeCuts1.removeDegenerates
                    (
                        10,              // max iterations
                        surf1,
                        querySurf2,
                        surf1PointTol,
                        points1         // work array
                    );

                if (nIters1 != 0)
                {
                    // Update geometric quantities
                    surf1.movePoints(points1);
                    hasMoved1 = true;
                }
            }
        }

        label nIters2 = 0;
        {
            triSurfaceSearch querySurf1(surf1);

            scalarField surf2PointTol
            (
                1e-3*edgeIntersections::minEdgeLength(surf2)
            );

            // Determine raw intersections
            edgeCuts2 = edgeIntersections
            (
                surf2,
                querySurf1,
                surf2PointTol
            );

            // Shuffle a bit to resolve degenerate edge-face hits
            {
                pointField points2(surf2.points());

                nIters2 =
                    edgeCuts2.removeDegenerates
                    (
                        10,              // max iterations
                        surf2,
                        querySurf1,
                        surf2PointTol,
                        points2         // work array
                    );

                if (nIters2 != 0)
                {
                    // Update geometric quantities
                    surf2.movePoints(points2);
                    hasMoved2 = true;
                }
            }
        }

        if (nIters1 == 0 && nIters2 == 0)
        {
            break;
        }
    }

    return hasMoved1 || hasMoved2;
}


void Foam::immersedStl::calcEdgeCuts
(
    triSurface& surf1,
    triSurface& surf2,
    const bool perturb,
    edgeIntersections& edge1Cuts,
    edgeIntersections& edge2Cuts
)
{
    if (perturb)
    {
        intersectSurfaces
        (
            surf1,
            edge1Cuts,
            surf2,
            edge2Cuts
        );
    }
    else
    {
        triSurfaceSearch querySurf2(surf2);

        edge1Cuts = edgeIntersections
        (
            surf1,
            querySurf2,
            1e-3*edgeIntersections::minEdgeLength(surf1)
        );

        triSurfaceSearch querySurf1(surf1);

        edge2Cuts = edgeIntersections
        (
            surf2,
            querySurf1,
            1e-3*edgeIntersections::minEdgeLength(surf2)
        );
    }
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::immersedStl::immersedStl
(
    const polyMesh& pMesh,
    const immersedBoundaryObject& ibo,
    const dictionary& dict
)
:
    immersedShape(pMesh, ibo, dict),
    dict_(dict),
    fileName_(dict.lookup("file")),
    immersedMesh_
    (
        IOobject
        (
            fileName_,
            pMesh.time().constant(), // instance
            "triSurface", // local
            pMesh, // read from parent registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    )
{
    read(dict);

    List<face> faces;
    immersedMesh_.movePoints(immersedMesh_.points() - this->centreOfMass_);
    vector cor(dict.lookupOrDefault<vector>("centreOfRotation", Zero));
    if (ei_ != -1)
    {
        //- Move to the initial location of the
        immersedMesh_.movePoints
        (
            immersedMesh_.points() + cor
        );
        get2DPoints(faces); // extrusion is done here
        points0_ -= cor;

        immersedMesh_.movePoints(immersedMesh_.points() - cor);

        //- Correct orientation of the original mesh
        orientedSurface::orient
        (
            immersedMesh_,
            vector::zero,
            true
        );
    }
    else
    {
        refine3D(immersedMesh_);
    }
    immersedMesh_.cleanup(false);
    points0_ = immersedMesh_.points();
    tssPtr_.set
    (
        new triSurfaceSearch(refCast<const triSurface>(immersedMesh_))
    );
    patchPtr_.set
    (
        new standAlonePatch(immersedMesh_.faces(), immersedMesh_.points())
    );

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedStl::~immersedStl()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Make list of points axisymmetric and create face indexing
void Foam::immersedStl::get2DPoints(List<face>& faces)
{
    points0_.clear();
    faces.clear();
    const polyBoundaryMesh& bMesh = pMesh_.boundaryMesh();
    List<wordRe> names(bMesh.names());
    labelHashSet patches(bMesh.patchSet(names));

    triSurface bTriMesh(triSurfaceTools::triangulate(bMesh, patches));

    PtrList<triSurface> regionTriSurfaces(1);
    regionTriSurfaces.set(0, new triSurface(immersedMesh_));
    boolList onAxis(1, true);
//     createRegionTriSurfaces(immersedMesh_, regionTriSurfaces, onAxis);

    forAll(regionTriSurfaces, regionI)
    {
        triSurface& immersedObject(regionTriSurfaces[regionI]);
        triSurfaceSearch tss(immersedObject);

        vector centre = boundBox(immersedObject.points()).midpoint();

        //- Make sure faces are pointing towards the interior
        boolList inside = tss.calcInside(pointField(1, centre));
        orientedSurface::orient(immersedObject, centre, !inside[0]);

        edgeIntersections edge1Cuts;
        edgeIntersections edge2Cuts;

        // Calculate the points where the edges are cut by the other surface
        calcEdgeCuts
        (
            immersedObject,
            bTriMesh,
            false,
            edge1Cuts,
            edge2Cuts
        );

        // Determine intersection edges from the edge cuts
        surfaceIntersection inter
        (
            immersedObject,
            edge1Cuts,
            bTriMesh,
            edge2Cuts
        );

        pointField points(inter.cutPoints());

        if (!returnReduce(points.size(), sumOp<label>()))
        {
            FatalErrorInFunction
                << "No intersection points were found."
                << "The specified file is either fully inside the" << nl
                << "mesh or fully outside" << endl
                << abort(FatalError);
        }

        pointField oldPoints(points);
        {
            //- Remove duplicate points and flatten mesh
            label pi = 0;
            forAll(oldPoints, i)
            {
                if (oldPoints[i][ei_] > 0.0)
                {
                    points[pi] = oldPoints[i];
                    points[pi][ei_] = 0.0;
                    pi++;
                }
            }
            points.resize(pi);
        }

        if (Pstream::parRun())
        {
            List<pointField> sendData(Pstream::nProcs());
            sendData[Pstream::myProcNo()] = points;
            Pstream::gatherList(sendData);
            Pstream::scatterList(sendData);

            points.clear();
            forAll(sendData, i)
            {
                points.append(sendData[i]);
            }
        }

        // Sort the points
        sortPointsPolar(points);

        // Add axis points if axisymmetric
        if (onAxis[regionI])
        {
            addAxisPoints(points);
        }
        // Sort the points
        sortPointsPolar(points);

        // Refine/coarsen region points
        coarsenRefine2D(points);

        List<face> subFaces;
        if (ai_ != -1)
        {
            extrudeAxi(points, subFaces, onAxis[regionI]);
        }
        else
        {
            extrude2D(points, subFaces);
        }

        label startI = points0_.size();
        forAll(subFaces, facei)
        {
            subFaces[facei][2] += startI;
            subFaces[facei][1] += startI;
            subFaces[facei][0] += startI;
        }

        points0_.append(points);
        faces.append(subFaces);
    }
}


void Foam::immersedStl::createRegionTriSurfaces
(
    const triSurface& triMesh,
    PtrList<triSurface>& regionTriSurfaces,
    boolList& onAxis
) const
{
    onAxis.clear();
    regionTriSurfaces.clear();

    PackedBoolList visited(triMesh.size(), false);
    visited.set(0, true);

    label startFace = 0;
    bool allVisited = false;
    while (!allVisited)
    {
        bool onAxisi = false;
        label fi = 0;
        label pi = 0;

        PackedBoolList pointAdded(triMesh.points().size(), false);
        PackedBoolList faceAdded(visited.size(), false);

        labelList pointMap(triMesh.points().size(), -1);
        labelList faceMap(faceAdded.size());

        pointField points(pointMap.size());

        //- Set starting face
        {
            vector maxP = Zero;
            vector minP = Zero;

            faceAdded.set(startFace, true);
            visited.set(startFace, true);
            faceMap[fi++] = startFace;
            forAll(triMesh[startFace], pj)
            {
                label pointj = triMesh[startFace][pj];

                const vector& p = triMesh.points()[pointj];
                maxP = max(maxP, p);
                minP = min(minP, p);

                pointAdded.set(pointj, true);
                points[pi] = p;
                pointMap[pointj] = pi++;
            }

            bool allNegative = maxP[yi_]*minP[yi_] < 0;
            allNegative =
                allNegative || maxP[ei_]*minP[ei_] < 0;
            onAxis = allNegative;
        }

        // Iterate over surface until all connecting faces have been added
        bool expanding = true;
        while (expanding)
        {
            expanding = false;

            forAll(triMesh, facei)
            {
                if (faceAdded.get(facei))
                {
                    forAll(triMesh.faceFaces()[facei], fj)
                    {
                        label facej = triMesh.faceFaces()[facei][fj];
                        if (!faceAdded.get(facej))
                        {
                            expanding = true;

                            faceAdded.set(facej, true);
                            faceMap[fi++] = facej;

                            vector maxP = Zero;
                            vector minP = Zero;
                            forAll(triMesh[facej], pj)
                            {
                                label pointj = triMesh[facej][pj];
                                const vector p = triMesh.points()[pointj];
                                maxP = max(maxP, p);
                                minP = min(minP, p);

                                if (!pointAdded.get(pointj))
                                {
                                    pointAdded.set(pointj, true);
                                    points[pi] = p;
                                    pointMap[pointj] = pi++;
                                }
                            }
                            bool allNegative = maxP[yi_]*minP[yi_] < 0;
                            allNegative =
                                allNegative || maxP[ei_]*minP[ei_] < 0;
                            onAxisi = onAxisi || allNegative;
                            visited.set(facej, true);
                        }
                    }
                }
            }
        }
        faceMap.resize(fi);
        points.resize(pi);

        List<labelledTri> faces(faceMap.size(), labelledTri(0, 0, 0, 0));

        forAll(faceMap, facei)
        {
            label facej = faceMap[facei];
            forAll(triMesh[facej], pj)
            {
                label pointj = triMesh[facej][pj];
                faces[facei][pj] = pointMap[pointj];
            }
        }

        regionTriSurfaces.append(new triSurface(faces, points));
        onAxis.append(onAxisi);

        allVisited = true;
        forAll(visited, i)
        {
            if (!visited.get(i))
            {
                startFace = i;
                allVisited = false;
                break;
            }
        }
    }
    Info<< "Found " << regionTriSurfaces.size() << " regions" << endl;
}


void Foam::immersedStl::coarsenRefine2D(pointField& points)
{
    //- Add points to meet goal mesh size
    pointField newPoints;
    label newi = 0;
    label nx = ai_ == -1 ? points.size() : points.size() - 1;
    for (label i = 0; i < nx; i++)
    {
        label j = (i + 1) % points.size();

        vector dx = points[j] - points[i];
        scalar magDx(mag(dx));
        if (magDx/dx_ > 1)
        {
            label nNew = label(magDx/dx_) + 1;
            newPoints.resize(newPoints.size() + nNew);

            vector newDx = 1.0/scalar(nNew + 2)*dx;
            for (label I = 1; I <= nNew; I++)
            {
                newPoints[newi++] = points[i] + newDx*scalar(I + 1);
            }
        }
    }
    points.append(newPoints);
    sortPointsPolar(points);

    pointField origPoints(points);

    // Points to remove
    PackedBoolList pointsToRemove(points.size(), false);

    // Starting point
    label i = 0;
    label j = 1;

    //- Remove points to meet goal mesh size
    bool good = true;
    while (good)
    {
        scalar dx = mag(points[i] - points[j]);
        if (dx < dx_)
        {
            pointsToRemove.set(j++);
        }
        else
        {
            i = j;
            j = ai_ == -1 ? (i + 1) % points.size() : i + 1;
        }

        if (j == points.size() && ai_ != -1)
        {
            good = false;
        }
        else if ( i > 2 && j == 2) // Stop after looping
        {
            good = false;
        }
    }

    // Do not remove the axis points
    pointsToRemove.set(0, false);
    if (ai_ != -1)
    {
        pointsToRemove.set(points.size()-1, false);
    }

    points.resize(points.size() - pointsToRemove.count());

    i = 0;
    forAll(origPoints, pi)
    {
        if (!pointsToRemove.get(pi))
        {
            points[i++] = origPoints[pi];
        }
    }
}

void Foam::immersedStl::addAxisPoints(pointField& points) const
{
    // Add axis points if axisymmetric
    if (ri_ != -1)
    {
        pointField origPoints(points);

        //- Add points on the axis
        points.resize(points.size() + 2);

        label pi = 0;
        //- Extrapolate to first point on axis
        {
            vector diff = origPoints[1] - origPoints[0];
            scalar x0 = origPoints[0][xi_];
            scalar y0 = origPoints[0][yi_];
            point xyz(Zero);

            if (mag(diff[ri_]) > small)
            {
                xyz[ai_] = x0 - y0*diff[xi_]/stabilise(diff[yi_], small);
            }
            else
            {
                xyz[xi_] = x0;
            }
            if (mag(origPoints[0] - xyz) > small)
            {
                points[pi++] = xyz;
            }
        }

        forAll(origPoints, i)
        {
            points[pi++] = origPoints[i];
        }

        //- Extrapolate to last point on axis
        {
            label nxi = origPoints.size() - 1;
            vector diff = origPoints[nxi] - origPoints[nxi-1];
            scalar x0 = origPoints[nxi][xi_];
            scalar y0 = origPoints[nxi][yi_];
            point xyz(Zero);

            if (mag(diff[ri_]) > small)
            {
                xyz[xi_] = x0 - y0*diff[xi_]/stabilise(diff[yi_], small);
            }
            else
            {
                xyz[xi_] = x0;
            }
            if (mag(origPoints.last() - xyz) > small)
            {
                points[pi++] = xyz;
            }
        }
        points.resize(pi);
    }
}


Foam::labelList
Foam::immersedStl::calcInside(const pointField& points) const
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

    boolList inside;
    if (ei_ != -1)
    {
        inside = tssPtr_->calcInside(object_.inverseTransform(validPoints));
    }
    else
    {
        inside = tssPtr_->calcInside(validPoints);
    }

    forAll(inside, i)
    {
        if(!inside[i]) // Flipped orientation
        {
            insidePoints[pi++] = insidePoints[i];
        }
    }
    insidePoints.resize(pi);
    return insidePoints;
}



bool Foam::immersedStl::inside(const point& pt) const
{
    if (bb_.contains(pt))
    {
        pointField pts(1, pt);
        if (ei_ != -1)
        {
            return tssPtr_->calcInside(object_.inverseTransform(pts))[0];
        }
        else
        {
            return tssPtr_->calcInside(pts)[0];
        }
    }
    return false;
}
// ************************************************************************* //
