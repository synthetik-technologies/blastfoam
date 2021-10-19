/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    2D Polygon intersection algorithms

Author
    Martin Beaudoin, Hydro-Quebec, (2008)

\*---------------------------------------------------------------------------*/

#include "boundBox.H"
#include "plane.H"

// Point in polygon algorithm
#include "HormannAgathos.H"

// Polygon clipping algorithm
#include "SutherlandHodgman.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// Computation of the intersection area between 2 2D polygons.
//
// Two following algorithms are considered here:
//   Sutherland-Hodgman algorithm : a classic, but polygons need to be convex
//   Greiner-Hormann algorithm    : very nice algorithm; accept both convex
//                                  and concave polygons
//
//   We also rely here on a very efficient algorithm for the point in
//   polygon problem for arbitrary polygons. That algorithm was proposed
//   by Hormann and Agathos.

template<class MasterPatch, class SlavePatch>
scalar
newGGIInterpolation<MasterPatch, SlavePatch>::polygonIntersection
(
    const List<point2D>& poly1,
    const List<point2D>& poly2
) const
{
    // Using pointers because clipping and subject may be swapped
    // by the algorithm.  HJ, 24/Oct/2008

     // The first polygon will be the clipping polygon
    const List<point2D>* clippingPolygon = &poly1;

    // Yhe second polygon will be the subject polygon; the one that
    // is going to be clipped
    const List<point2D>* subjectPolygon = &poly2;

    // Empty list so we can detect weird cases later
    List<point2D> clippedPolygon;
    scalar intersectionArea = 0.0;

    // First, let's get rid of the obvious:
    //  1: Neither polygons intersect one another
    //     --> intersection area == 0.0 : that should not happen
    //  2: If both polygons completely overlap one another
    //     ---> subjectPolygon is the intersection polygon
    //  3: If clippingPolygon totally enclosed subjectPolygon
    //     ---> subjectPolygon is the intersecting polygon
    //  4: If subjectPolygon totally enclosed clippingPolygon
    //     --> clippingPolygon is the intersecting polygon
    //  5: Otherwise, we have partial or full intersection
    //
    //  For this, we first detect if the vertices of the subject polygon
    //  are inside or outside of the clipping polygon
    //  This is a quick intersection test...

    // Keep track of who is inside or outside
    List<bool> subjectVertexInside(subjectPolygon->size());

    insideOutside statusInOut =
        isVertexInsidePolygon
        (
            *clippingPolygon,
            *subjectPolygon,
            subjectVertexInside
        );

    // We check right away if statusInOut == ALL_OUTSIDE
    // Sometimes, it is just that the clipping polygon is inside the
    // subject polygon instead So let's explore this situation, just
    // in case
    if (statusInOut == ALL_OUTSIDE)
    {
        // We need to check if subject is not completely or partially
        // enclosing clipping instead

        clippingPolygon = &poly2;
        subjectPolygon  = &poly1;

        subjectVertexInside.setSize(subjectPolygon->size());
        statusInOut =
            isVertexInsidePolygon
            (
                *clippingPolygon,
                *subjectPolygon,
                subjectVertexInside
            );

    }

    switch(statusInOut)
    {
        case ALL_INSIDE:
        {
            clippedPolygon = *subjectPolygon;
            break;
        }
        case ALL_OUTSIDE:
        case PARTIALLY_OVERLAPPING:
        default:
        {
            // Compute the intersection

            // If by any chance, we have reached a situation where the
            // intersection is really zero, it is because the quick
            // reject tests have missed something. The intersection
            // area will be 0.0, and the calling function should at
            // least check for this, and report loudly...
            clippedPolygon =
                clipPolygon2DSutherlandHodgman
                (
                    *clippingPolygon,
                    *subjectPolygon
                );// For convex polygons

#if 0
            // This is the next candidate to code if the Sutherland
            // Hodgman algorithm encounters concave polygons...
            clippedPolygon =
                clipPolygon2DGreinerHormann
                (
                    *clippingPolygon,
                    *subjectPolygon,
                    subjectVertexInside
                );   // For arbitrary polygons
#endif
            break;
        }
    }

    // Compute the area of clippedPolygon if we indeed do have a
    // clipped polygon; otherwise, the computed area stays at 0.0;
    if (clippedPolygon.size() > 2)
    {
        // We are only interested in the absolute value
        intersectionArea = mag(area2D(clippedPolygon));
    }

    if (debug)
    {
        // Check against tolerances
        scalar clippingArea  = area2D(*clippingPolygon);
        scalar subjectArea   = area2D(*subjectPolygon);

        if
        (
            mag(intersectionArea/clippingArea) < areaErrorTol_()
         || mag(intersectionArea/subjectArea)  < areaErrorTol_()
        )
        {
            WarningIn
            (
                "newGGIInterpolation<MasterPatch, SlavePatch>::"
                "polygonIntersection"
            )   << "Intersection might be wrong.  Clipping side "
                << intersectionArea/clippingArea << " subject: "
                << intersectionArea/subjectArea << endl;
        }
    }

    return intersectionArea;
}


// Compute the point in polygon problem using the computation of the
// winding number. This technique is based on a paper by Hormann and
// Agathos.
//
// This is based on the winding number technique, but optimized in
// order to only evaluate quarter revolutions instead of the whole
// arcos/sqrt basic algorithm.  When the GGI weighting factors will
// have to be recomputed often for moving meshes, this performance
// will be useful.  We can also go back to the classical winding
// number algorithm if need be Notice: The list subjectVertexInside
// will return a boolean marking if a point from the subject polygon
// is inside or outside the clipping polygon. A point is "inside"
// (value == true) if it is inside the clipping polygon, or on a
// vertex or edge.  This information will be very useful for the
// Sutherland Hodgman algo.
template<class MasterPatch, class SlavePatch>
typename Foam::newGGIInterpolation<MasterPatch, SlavePatch>::insideOutside
newGGIInterpolation<MasterPatch, SlavePatch>::isVertexInsidePolygon
(
    const List<point2D>& clippingPolygon,
    const List<point2D>& subjectPolygon,
    List<bool>& subjectVertexInside
) const
{
    insideOutside retValue = ALL_OUTSIDE;

    // The class HormannAgathos implements the algorithm
    // We use distErrorTol_ to evaluate a distance factor called
    // epsilon.  That epsilon factor will be used to detect if a point
    // is on a vertex or an edge.
    const scalar distErrorTol = sqrt(areaErrorTol_());

    HormannAgathos pip(clippingPolygon, distErrorTol);

    // Counter
    label nbrsOutside = 0;

    // Iterate over all subject vertices, determining if they are
    // inside/outside the clipping polygon
    forAll(subjectPolygon, sPI)
    {
        switch (pip.evaluate(subjectPolygon[sPI]))
        {
            case HormannAgathos::POINT_OUTSIDE:
                nbrsOutside++;
                subjectVertexInside[sPI] = false;
                break;
            case HormannAgathos::POINT_ON_VERTEX:
            case HormannAgathos::POINT_ON_EDGE:
            case HormannAgathos::POINT_INSIDE:
            default:
                // Vertex, Edge or Inside: We are inside!
                subjectVertexInside[sPI] = true;
                break;
        }
    }

    // Let's do the inventory now...
    if (nbrsOutside == 0)
    {
        retValue = ALL_INSIDE;
    }
    else if (nbrsOutside < subjectPolygon.size())
    {
        retValue = PARTIALLY_OVERLAPPING;
    }
    // else, all the points are outside, which is not necessarily a
    // problem if the subject Polygon is enclosing partially or
    // completely the clipping polygon instead

    return retValue;
}


// This is an implementation of the Sutherland Hodgman algorithm:
// Reentrant Polygon Clipping, Sutherland, Hodgman, Communications of
// the ACM, 1974
//
// Wikipedia has a very simple Pseudo-Code example of this rather
// simple algorithm http://en.wikipedia.org/wiki/Sutherland-Hodgeman.
//
// The subject polygon will be clipped by the clipping polygon. The
// list of boolean values subjectVertexInside provide the bonus
// information if a subject vertex is either "inside" or "outside"
// the clipping polygon. A vertex is "inside" if it is inside, on a
// vertex, or on an edge. Otherwise, the vertex is "outside"
template<class MasterPatch, class SlavePatch>
List<point2D>
newGGIInterpolation<MasterPatch, SlavePatch>::clipPolygon2DSutherlandHodgman
(
    const List<point2D>& clippingPolygon,
    const List<point2D>& subjectPolygon
) const
{
    // Create Sutherland-Hodgman intersector, with a distance tolerance
    // factor for intersection computation
    // We use a tolerance that is consistant with  the quick reject tests
    return SutherlandHodgman
    (
        clippingPolygon,
        subjectPolygon,
        sqrt(areaErrorTol_()) // = distErrorTol
    ).evaluate();
}


// Compute the area of a polygon
// See:  http://mathworld.wolfram.com/PolygonArea.html
template<class MasterPatch, class SlavePatch>
scalar newGGIInterpolation<MasterPatch, SlavePatch>::area2D
(
    const List<point2D>& polygon
) const
{
    // For a non-self-intersecting (simple) polygon with n vertices,
    // the area is :
    // A = 0.5 * sum(x(i)y(i+1) - x(i+1)y(i)  , i=1 to n; where n+1 = 0;
    scalar area = 0;

    // We start with last term from Wolfram equation
    label indexI   = polygon.size() - 1;
    label indexIp1 = indexI;

    for (label i = 0; i < polygon.size(); i++)
    {
        indexI = indexIp1;
        indexIp1 = i;

        area += polygon[indexI].x()*polygon[indexIp1].y()
              - polygon[indexIp1].x()*polygon[indexI].y();
    }

    // NB: THe area can be positive or negative
    return area/2.0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
