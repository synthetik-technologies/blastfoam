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

Author
    Martin Beaudoin, Hydro-Quebec, (2008)

\*---------------------------------------------------------------------------*/

#include "SutherlandHodgman.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Compute the intersection between two line segments
// See:
//     http://mathworld.wolfram.com/Line-LineIntersection.html
//     http://local.wasp.uwa.edu.au/~pbourke/geometry/lineline2d/
//
//  We are only interested in the actual interpolation where ua and ub are
//  both in the range [0, 1]. If outside of this range, we simply return false.
bool Foam::SutherlandHodgman::lineSegmentIntersection
(
    const point2D& p1,  // Segment 1 [p1, p2]
    const point2D& p2,
    const point2D& p3,  // Segment 2 [p3, p4]
    const point2D& p4,
    point2D& intersectionPoint
) const
{
    bool segmentDoIntersect = false;

    intersectionPoint[0] = 0.0;
    intersectionPoint[1] = 0.0;

    vector2D delta21(p2 - p1);
    vector2D delta43(p4 - p3);
    vector2D delta13(p1 - p3);

    scalar denom = delta43[1]*delta21[0] - delta43[0]*delta21[1];

    // Comparison with 0 and 1...  Floating point precision
    // issue. Cannot compare with 0 and 1 exactly.  VSMALL is too
    // strict... we need to take into account that we tolerate an
    // epsilon error all along so in the worst case, if a point P is
    // separated by only epsilon from an edge, the smallest
    // determinant should be of the order of 2*(2*epsilon*2epsilon) =
    // 8*epsilon^2
    //
    // Epsilon is adjusted to the smallest edge length between [p1,p2]
    // and [p3,p4]

    // Accounting for double precision tolerance.  HJ, 16/Jan/2009
    scalar epsilon =
        Foam::max
        (
            Foam::min(mag(delta21), mag(delta43))*intersectSegDistTol_,
            SMALL
        );

    scalar tolFactor = 8.0*sqr(epsilon);

    if (mag(denom) > tolFactor)
    {
        scalar ua = (delta43[0]*delta13[1] - delta43[1]*delta13[0])/denom;
        scalar ub = (delta21[0]*delta13[1] - delta21[1]*delta13[0])/denom;

        //- Adjust tolFactor because of division:
        //  (y +- epsilon)/(x +- epsilon) = y/x +- 2*epsilon (rough estimation)
        tolFactor *= 2.0;

        // Check for end points outside of range 0..1
        //if (ua >= 0.0 && ua <= 1.0 || ub >= 0.0 && ub <= 1.0)
        // One of the intersection must be inside of one of the segment
        // We don't want both of them at the same time, because each edge
        // should be considered as infinite line here. But still, the
        // intersection must be within the interval [0,1] for one of
        // the edge
        if (ua >= -tolFactor && ua <= (1.0 + tolFactor))
        {
            // This is the intersection we want
            intersectionPoint = p1 + ua*delta21;
            segmentDoIntersect = true;

            if (mag(intersectionPoint - (p3 + ub*delta43)) > epsilon)
            {
                WarningIn
                (
                    "Foam::SutherlandHodgman::lineSegmentIntersection()"
                )   << "ua does not match with ub: delta: "
                    << mag(intersectionPoint - (p3 + ub*delta43))
                    << " : epsilon: " << epsilon
                    << endl;
            }
        }
        else if (ub >= -tolFactor && ub <= (1.0 + tolFactor))
        {
            // This is the intersection we want
            intersectionPoint = p3 + ub*delta43;
            segmentDoIntersect = true;

            if (mag(intersectionPoint - (p1 + ua*delta21)) > epsilon)
            {
                WarningIn
                (
                    "Foam::SutherlandHodgman::lineSegmentIntersection()"
                )   << "ub does not match with ua: delta: "
                    << mag(intersectionPoint - (p1 + ua*delta21))
                    << " : epsilon: " << epsilon
                    << endl;
            }
        }
    }

    // Parallel lines. We don't go any further
    return segmentDoIntersect;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SutherlandHodgman::SutherlandHodgman
(
    const List<point2D>& clippingPolygon,
    const List<point2D>& subjectPolygon,
    const scalar& intersectSegDistTol
)
:
    subjectPolygon_(subjectPolygon),
    clippingPolygon_(clippingPolygon),
    currentClipEdgeP1_(clippingPolygon.size() - 1),
    currentClipEdgeP2_(0),
    intersectSegDistTol_(intersectSegDistTol)
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::List<Foam::point2D> Foam::SutherlandHodgman::evaluate()
{
    DynamicList<point2D, 8> clippedPolygon(0);
    DynamicList<bool, 8> clippedVertexInside(0);
    point2D interSectPt;

    label S = subjectPolygon_.size() - 1;

    forAll (subjectPolygon_, sPI)
    {
        label P  = sPI;

        // if P is inside or visible from cutting edge
        if (isPointVisibleFromCuttingSegment(subjectPolygon_[P]))
        {
            // if S is inside or visible from cutting edge
            if (isPointVisibleFromCuttingSegment(subjectPolygon_[S]))
            {
                clippedPolygon.append(subjectPolygon_[P]);  // Output P
            }
            else
            {
                // Output intersection between S and P
                if
                (
                    clipSegment
                    (
                        subjectPolygon_[S],
                        subjectPolygon_[P],
                        interSectPt
                    )
                )
                {
                    // We do intersect this edge

                    // Output intersection point
                    clippedPolygon.append(interSectPt);

                    // Build the inside/outside list
                    clippedVertexInside.append(true);
                }
                else
                {
//                     WarningIn("List<point2D> SutherlandHodgman::evaluate()")
//                         << "missed intersection... recheck algorithm, 1"
//                         << endl;
                }

                // Output P
                clippedPolygon.append(subjectPolygon_[P]);
            }
        }
        // if S is inside or visible from cutting edge
        else if (isPointVisibleFromCuttingSegment(subjectPolygon_[S]))
        {
            if
            (
                clipSegment
                (
                    subjectPolygon_[P],
                    subjectPolygon_[S],
                    interSectPt
                )
            )
            {
                // We do intersect this edge

                // Output intersection point
                clippedPolygon.append(interSectPt);

                // Build the inside/outside list
                clippedVertexInside.append(true);
            }
            else
            {
//                 WarningIn("List<point2D> SutherlandHodgman::evaluate()")
//                     << "missed intersection... recheck algorithm, 2"
//                     << endl;
            }
        }

        S = P;
    }

    // Next clipping edge
    currentClipEdgeP1_ = currentClipEdgeP2_;
    currentClipEdgeP2_++;

    // Update, clean up and return
    subjectPolygon_.transfer(clippedPolygon.shrink());

    // Are we done?
    if (currentClipEdgeP2_ < clippingPolygon_.size())
    {
        // Clip against next edge, re-entrant.
        evaluate();
    }

    // We are done
    return subjectPolygon_; // Which is now totally clipped
}


// ************************************************************************* //
