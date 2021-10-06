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

#include "octreeDataFace.H"
#include "labelList.H"
#include "polyMesh.H"
#include "octree.H"
#include "polyPatch.H"
#include "triangleFuncs.H"
#include "linePointRef.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::octreeDataFace, 0);

Foam::scalar Foam::octreeDataFace::tol = 1e-6;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::octreeDataFace::calcBb()
{
    allBb_.setSize(meshFaces_.size());
    allBb_ = treeBoundBox::invertedBox;

    forAll (meshFaces_, i)
    {
        // Update bb of face
        treeBoundBox& myBb = allBb_[i];

        const face& f = mesh_.faces()[meshFaces_[i]];

        forAll(f, faceVertexI)
        {
            const point& coord = mesh_.points()[f[faceVertexI]];

            myBb.min() = min(myBb.min(), coord);
            myBb.max() = max(myBb.max(), coord);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from selected mesh faces
Foam::octreeDataFace::octreeDataFace
(
    const primitiveMesh& mesh,
    const labelList& meshFaces,
    const treeBoundBoxList& allBb
)
:
    mesh_(mesh),
    meshFaces_(meshFaces),
    allBb_(allBb)
{}


// Construct from selected mesh faces. Bounding box calculated.
Foam::octreeDataFace::octreeDataFace
(
    const primitiveMesh& mesh,
    const labelList& meshFaces
)
:
    mesh_(mesh),
    meshFaces_(meshFaces),
    allBb_(meshFaces_.size())
{
    // Generate tight fitting bounding box
    calcBb();
}


// Construct from selected mesh faces
Foam::octreeDataFace::octreeDataFace
(
    const primitiveMesh& mesh,
    const UList<const labelList*>& meshFaceListPtrs,
    const UList<const treeBoundBoxList*>& bbListPtrs
)
:
    mesh_(mesh),
    meshFaces_(0),
    allBb_(0)
{
    label faceI = 0;

    forAll(meshFaceListPtrs, listI)
    {
        faceI += meshFaceListPtrs[listI]->size();
    }

    meshFaces_.setSize(faceI);
    allBb_.setSize(faceI);

    faceI = 0;

    forAll(meshFaceListPtrs, listI)
    {
        const labelList& meshFaces = *meshFaceListPtrs[listI];
        const treeBoundBoxList& allBb = *bbListPtrs[listI];

        forAll(meshFaces, meshFaceI)
        {
            meshFaces_[faceI] = meshFaces[meshFaceI];
            allBb_[faceI] = allBb[meshFaceI];
            faceI++;
        }
    }
}


// Construct from selected mesh faces. Bounding box calculated.
Foam::octreeDataFace::octreeDataFace
(
    const primitiveMesh& mesh,
    const UList<const labelList*>& meshFaceListPtrs
)
:
    mesh_(mesh),
    meshFaces_(0)
{
    label faceI = 0;

    forAll(meshFaceListPtrs, listI)
    {
        faceI += meshFaceListPtrs[listI]->size();
    }

    meshFaces_.setSize(faceI);

    faceI = 0;

    forAll(meshFaceListPtrs, listI)
    {
        const labelList& meshFaces = *meshFaceListPtrs[listI];

        forAll(meshFaces, meshFaceI)
        {
            meshFaces_[faceI++] = meshFaces[meshFaceI];
        }
    }

    // Generate tight fitting bounding box
    calcBb();
}


// Construct from all faces in polyPatch. Bounding box calculated.
Foam::octreeDataFace::octreeDataFace(const polyPatch& patch)
:
    mesh_(patch.boundaryMesh().mesh()),
    meshFaces_(patch.size())
{
    forAll(patch, patchFaceI)
    {
        meshFaces_[patchFaceI] = patch.start() + patchFaceI;
    }

    // Generate tight fitting bounding box
    calcBb();
}


// Construct from primitiveMesh. Inserts all boundary faces.
Foam::octreeDataFace::octreeDataFace(const primitiveMesh& mesh)
:
    mesh_(mesh),
    meshFaces_(0),
    allBb_(0)
{
    // Size storage
    meshFaces_.setSize(mesh_.nFaces() - mesh_.nInternalFaces());

    // Set info for all boundary faces.
    label boundaryFaceI = 0;

    for (label faceI = mesh_.nInternalFaces(); faceI < mesh_.nFaces(); faceI++)
    {
        meshFaces_[boundaryFaceI++] = faceI;
    }

    // Generate tight fitting bounding box
    calcBb();
}


// Construct as copy
Foam::octreeDataFace::octreeDataFace(const octreeDataFace& shapes)
:
    mesh_(shapes.mesh()),
    meshFaces_(shapes.meshFaces()),
    allBb_(shapes.allBb())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::octreeDataFace::~octreeDataFace()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::octreeDataFace::getSampleType
(
    const octree<octreeDataFace>& oc,
    const point& sample
) const
{
    // Need to determine whether sample is 'inside' or 'outside'
    // Done by finding nearest face. This gives back a face which is
    // guaranteed to contain nearest point. This point can be
    // - in interior of face: compare to face normal
    // - on edge of face: compare to edge normal
    // - on point of face: compare to point normal
    // Unfortunately the octree does not give us back the intersection point
    // or where on the face it has hit so we have to recreate all that
    // information.

    treeBoundBox tightest(treeBoundBox::greatBox);
    scalar tightestDist(treeBoundBox::great);
    // Find nearest face to sample
    label index = oc.findNearest(sample, tightest, tightestDist);

    if (index == -1)
    {
        FatalErrorIn
        (
            "octreeDataFace::getSampleType"
            "(octree<octreeDataFace>&, const point&)"
        )   << "Could not find " << sample << " in octree."
            << abort(FatalError);
    }


    // Get actual intersection point on face
    label faceI = meshFaces_[index];

    if (debug & 2)
    {
        Pout<< "getSampleType : sample:" << sample
            << " nearest face:" << faceI;
    }

    const face& f = mesh_.faces()[faceI];

    const pointField& points = mesh_.points();

    pointHit curHit = f.nearestPoint(sample, points);

    //
    // 1] Check whether sample is above face
    //

    if (curHit.hit())
    {
        // Simple case. Compare to face normal.

        if (debug & 2)
        {
            Pout<< " -> face hit:" << curHit.hitPoint()
                << " comparing to face normal " << mesh_.faceAreas()[faceI]
                << endl;
        }
        return octree<octreeDataFace>::getVolType
        (
            mesh_.faceAreas()[faceI],
            sample - curHit.hitPoint()
        );
    }

    if (debug & 2)
    {
        Pout<< " -> face miss:" << curHit.missPoint();
    }

    //
    // 2] Check whether intersection is on one of the face vertices or
    //    face centre
    //

    scalar typDim = sqrt(mag(mesh_.faceAreas()[faceI])) + VSMALL;

    forAll(f, fp)
    {
        if ((mag(points[f[fp]] - curHit.missPoint())/typDim) < tol)
        {
            // Face intersection point equals face vertex fp

            // Calculate point normal (wrong: uses face normals instead of
            // triangle normals)
            const labelList& myFaces = mesh_.pointFaces()[f[fp]];

            vector pointNormal(vector::zero);

            forAll(myFaces, myFaceI)
            {
                if (myFaces[myFaceI] >= mesh_.nInternalFaces())
                {
                    vector n = mesh_.faceAreas()[myFaces[myFaceI]];
                    n /= mag(n) + VSMALL;

                    pointNormal += n;
                }
            }

            if (debug & 2)
            {
                    Pout<< " -> face point hit :" << points[f[fp]]
                        << " point normal:" << pointNormal
                        << " distance:"
                        << mag(points[f[fp]] - curHit.missPoint())/typDim
                        << endl;
            }
            return octree<octreeDataFace>::getVolType
            (
                pointNormal,
                sample - curHit.missPoint()
            );
        }
    }
    if ((mag(mesh_.faceCentres()[faceI] - curHit.missPoint())/typDim) < tol)
    {
        // Face intersection point equals face centre. Normal at face centre
        // is already average of face normals

        if (debug & 2)
        {
            Pout<< " -> centre hit:" << mesh_.faceCentres()[faceI]
                << " distance:"
                << mag(mesh_.faceCentres()[faceI] - curHit.missPoint())/typDim
                << endl;
        }

        return octree<octreeDataFace>::getVolType
        (
            mesh_.faceAreas()[faceI],
            sample - curHit.missPoint()
        );
    }


    //
    // 3] Get the 'real' edge the face intersection is on
    //

    const labelList& myEdges = mesh_.faceEdges()[faceI];

    forAll(myEdges, myEdgeI)
    {
        const edge& e = mesh_.edges()[myEdges[myEdgeI]];

        pointHit edgeHit = line<point, const point&>
        (
            points[e.start()],
            points[e.end()]
        ).nearestDist(sample);


        if ((mag(edgeHit.rawPoint() - curHit.missPoint())/typDim) < tol)
        {
            // Face intersection point lies on edge e

            // Calculate edge normal (wrong: uses face normals instead of
            // triangle normals)
            const labelList& myFaces = mesh_.edgeFaces()[myEdges[myEdgeI]];

            vector edgeNormal(vector::zero);

            forAll(myFaces, myFaceI)
            {
                if (myFaces[myFaceI] >= mesh_.nInternalFaces())
                {
                    vector n = mesh_.faceAreas()[myFaces[myFaceI]];
                    n /= mag(n) + VSMALL;

                    edgeNormal += n;
                }
            }

            if (debug & 2)
            {
                Pout<< " -> real edge hit point:" << edgeHit.rawPoint()
                    << " comparing to edge normal:" << edgeNormal
                    << endl;
            }

            // Found face intersection point on this edge. Compare to edge
            // normal
            return octree<octreeDataFace>::getVolType
            (
                edgeNormal,
                sample - curHit.missPoint()
            );
        }
    }


    //
    // 4] Get the internal edge the face intersection is on
    //

    forAll(f, fp)
    {
        pointHit edgeHit =
            line<point, const point&>
            (
                points[f[fp]],
                mesh_.faceCentres()[faceI]
            ).nearestDist(sample);

        if ((mag(edgeHit.rawPoint() - curHit.missPoint())/typDim) < tol)
        {
            // Face intersection point lies on edge between two face triangles

            // Calculate edge normal as average of the two triangle normals
            const label fpPrev = f.rcIndex(fp);
            const label fpNext = f.fcIndex(fp);

            vector e = points[f[fp]] - mesh_.faceCentres()[faceI];
            vector ePrev = points[f[fpPrev]] - mesh_.faceCentres()[faceI];
            vector eNext = points[f[fpNext]] - mesh_.faceCentres()[faceI];

            vector nLeft = ePrev ^ e;
            nLeft /= mag(nLeft) + VSMALL;

            vector nRight = e ^ eNext;
            nRight /= mag(nRight) + VSMALL;

            if (debug & 2)
            {
                Pout<< " -> internal edge hit point:" << edgeHit.rawPoint()
                    << " comparing to edge normal "
                    << 0.5*(nLeft + nRight)
                    << endl;
            }

            // Found face intersection point on this edge. Compare to edge
            // normal
            return octree<octreeDataFace>::getVolType
            (
                0.5*(nLeft + nRight),
                sample - curHit.missPoint()
            );
        }
    }

    if (debug & 2)
    {
        Pout<< "Did not find sample " << sample
            << " anywhere related to nearest face " << faceI << endl
            << "Face:";

        forAll(f, fp)
        {
            Pout<< "    vertex:" << f[fp] << "  coord:" << points[f[fp]]
                << endl;
        }
    }

    // Can't determine status of sample with respect to nearest face.
    // Either
    // - tolerances are wrong. (if e.g. face has zero area)
    // - or (more likely) surface is not closed.

    return octree<octreeDataFace>::UNKNOWN;
}


bool Foam::octreeDataFace::overlaps
(
    const label index,
    const treeBoundBox& sampleBb
) const
{
    //return sampleBb.overlaps(allBb_[index]);

    //- Exact test of face intersecting bb

    // 1. Quick rejection: bb does not intersect face bb at all
    if (!sampleBb.overlaps(allBb_[index]))
    {
        return false;
    }

    // 2. Check if one or more face points inside
    label faceI = meshFaces_[index];

    const face& f = mesh_.faces()[faceI];

    const pointField& points = mesh_.points();

    forAll(f, fp)
    {
        if (sampleBb.contains(points[f[fp]]))
        {
            return true;
        }
    }

    // 3. Difficult case: all points are outside but connecting edges might
    // go through cube. Use triangle-bounding box intersection.
    const point& fc = mesh_.faceCentres()[faceI];

    forAll(f, fp)
    {
        const label fp1 = f.fcIndex(fp);

        bool triIntersects = triangleFuncs::intersectBb
        (
            points[f[fp]],
            points[f[fp1]],
            fc,
            sampleBb
        );

        if (triIntersects)
        {
            return true;
        }
    }

    return false;
}


bool Foam::octreeDataFace::contains(const label, const point&) const
{
    notImplemented
    (
        "octreeDataFace::contains(const label, const point&)"
    );

    return false;
}


bool Foam::octreeDataFace::intersects
(
    const label index,
    const point& start,
    const point& end,
    point& intersectionPoint
) const
{
    label faceI = meshFaces_[index];

    const face& f = mesh_.faces()[faceI];

    const vector dir(end - start);

    // Disable picking up intersections behind us.
    scalar oldTol = intersection::setPlanarTol(0.0);

    pointHit inter = f.ray
    (
        start,
        dir,
        mesh_.points(),
        intersection::algorithm::halfRay,
        intersection::direction::vector
    );

    intersection::setPlanarTol(oldTol);

    if (inter.hit() && inter.distance() <= mag(dir))
    {
        intersectionPoint = inter.hitPoint();

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::octreeDataFace::findTightest
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
Foam::scalar Foam::octreeDataFace::calcSign
(
    const label index,
    const point& sample,
    point& n
) const
{
    label faceI = meshFaces_[index];

    n = mesh_.faceAreas()[faceI];

    n /= mag(n) + VSMALL;

    vector vec = sample - mesh_.faceCentres()[faceI];

    vec /= mag(vec) + VSMALL;

    return n & vec;
}


// Calculate nearest point on/in shapei
Foam::scalar Foam::octreeDataFace::calcNearest
(
    const label index,
    const point& sample,
    point& nearest
) const
{
    label faceI = meshFaces_[index];

    const face& f = mesh_.faces()[faceI];

    pointHit nearHit = f.nearestPoint(sample, mesh_.points());

    nearest = nearHit.rawPoint();

    if (debug & 1)
    {
        const point& ctr = mesh_.faceCentres()[faceI];

        scalar sign = mesh_.faceAreas()[faceI] & (sample - nearest);

        Pout<< "octreeDataFace::calcNearest : "
            << "sample:" << sample
            << "  index:" << index
            << "  faceI:" << faceI
            << "  ctr:" << ctr
            << "  sign:" << sign
            << "  nearest point:" << nearest
            << "  distance to face:" << nearHit.distance()
            << endl;
    }
    return nearHit.distance();
}


// Calculate nearest point on/in shapei
Foam::scalar Foam::octreeDataFace::calcNearest
(
    const label index,
    const linePointRef& ln,
    point& linePt,
    point& shapePt
) const
{
    notImplemented
    (
        "octreeDataFace::calcNearest(const label, const linePointRef&"
        ", point&, point&)"
    );
    return GREAT;
}


void Foam::octreeDataFace::write(Ostream& os, const label index) const
{
    os << meshFaces_[index] << " " << allBb_[index];
}


// ************************************************************************* //
