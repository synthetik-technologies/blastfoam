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

#include "searchableAnnulusCylinder.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(searchableAnnulusCylinder, 0);
    addToRunTimeSelectionTable
    (
        searchableSurface,
        searchableAnnulusCylinder,
        dict
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::searchableAnnulusCylinder::coordinates() const
{
    tmp<pointField> tCtrs(new pointField(1, 0.5*(point1_ + point2_)));

    return tCtrs;
}


void Foam::searchableAnnulusCylinder::boundingSpheres
(
    pointField& centres,
    scalarField& radiusSqr
) const
{
    centres.setSize(1);
    centres[0] = 0.5*(point1_ + point2_);

    radiusSqr.setSize(1);
    radiusSqr[0] = Foam::magSqr(point1_-centres[0]) + Foam::sqr(outerRadius_);

    // Add a bit to make sure all points are tested inside
    radiusSqr += Foam::sqr(small);
}


Foam::tmp<Foam::pointField> Foam::searchableAnnulusCylinder::points() const
{
    tmp<pointField> tPts(new pointField(2));
    pointField& pts = tPts.ref();

    pts[0] = point1_;
    pts[1] = point2_;

    return tPts;
}


Foam::pointIndexHit Foam::searchableAnnulusCylinder::findNearest
(
    const point& sample,
    const scalar nearestDistSqr
) const
{
    pointIndexHit info(false, sample, -1);

    vector v(sample - point1_);

    // Decompose sample-point1 into normal and parallel component
    scalar parallel = (v & unitDir_);

    // Remove the parallel component and normalise
    v -= parallel*unitDir_;
    scalar magV = mag(v);

    if (magV < rootVSmall)
    {
        v = Zero;
    }
    else
    {
        v /= magV;
    }

    scalar distSqr = great;

    if (parallel <= 0)
    {
        // nearest is at point1 end cap. Clip to radius.
        scalar r = magV;
        if (magV < innerRadius_)
        {
            r = innerRadius_;
        }
        else if (magV > outerRadius_)
        {
            r = outerRadius_;
        }
        info.setPoint(point1_ + r*v);
        distSqr = magSqr(sample - info.rawPoint());

    }
    else if (parallel >= magDir_)
    {
        // nearest is at point2 end cap. Clip to radius.
        scalar r = magV;
        if (magV < innerRadius_)
        {
            r = innerRadius_;
        }
        else if (magV > outerRadius_)
        {
            r = outerRadius_;
        }
        info.setPoint(point2_ + r*v);
        distSqr = magSqr(sample - info.rawPoint());
    }
    else
    {
        // in between endcaps. Might either be nearer endcaps or cylinder wall

        // distance to endpoint: parallel or parallel-magDir
        // distance to cylinder wall: magV-radius

        // Nearest cylinder point
        if (magV < rootVSmall)
        {
            // Point exactly on centre line. Take any point on wall.
            vector e1 = point(1,0,0) ^ unitDir_;
            scalar magE1 = mag(e1);
            if (magE1 < small)
            {
                e1 = point(0,1,0) ^ unitDir_;
                magE1 = mag(e1);
            }
            e1 /= magE1;
            info.setPoint(sample + innerRadius_*e1);
            distSqr = magSqr(sample - info.rawPoint());
        }
        else if (magV <= innerRadius_)
        {
            info.setPoint(sample + (innerRadius_ - magV)*v);
            distSqr = magSqr(magV - innerRadius_);
        }
        else
        {
            const point cylPtI = sample + (innerRadius_ - magV)*v;
            const point cylPtO = sample + (outerRadius_ - magV)*v;

            // Project onto nearest endcap
            const point endPt =
                parallel < 0.5*magDir_
              ? point1_ + min(max(magV, innerRadius_), outerRadius_)*v
              : point2_ + min(max(magV, innerRadius_), outerRadius_)*v;

            const scalar endDist = magSqr(sample - endPt);
            const scalar cylDistI = magSqr(innerRadius_ - magV);
            const scalar cylDistO = magSqr(outerRadius_ - magV);
            if (endDist < cylDistI && endDist < cylDistO)
            {
                info.setPoint(endPt);
                distSqr = endDist;
            }
            else if (cylDistI < cylDistO)
            {
                info.setPoint(cylPtI);
                distSqr = cylDistI;
            }
            else
            {
                info.setPoint(cylPtO);
                distSqr = cylDistO;
            }
        }
    }

    if (distSqr < nearestDistSqr)
    {
        info.setHit();
        info.setIndex(0);
    }

    return info;
}


Foam::scalar Foam::searchableAnnulusCylinder::radius2(const point& pt) const
{
    const vector x = (pt - point1_) ^ unitDir_;
    return x & x;
}


// From http://www.gamedev.net/community/forums/topic.asp?topic_id=467789 -
// intersection of cylinder with ray
Foam::label Foam::searchableAnnulusCylinder::findLineAll
(
    const point& start,
    const point& end,
    pointIndexHit& nearO,
    pointIndexHit& nearI,
    pointIndexHit& farI,
    pointIndexHit& farO
) const
{
    nearO.setMiss();
    nearI.setMiss();
    farI.setMiss();
    farO.setMiss();

    vector point1Start(start - point1_);
    vector point2Start(start - point2_);
    vector point1End(end - point1_);

    // Quick rejection of complete vector outside endcaps
    scalar s1 = point1Start & unitDir_;
    scalar s2 = point1End & unitDir_;

    if ((s1 < 0 && s2 < 0) || (s1 > magDir_ && s2 > magDir_))
    {
        return 0;
    }

    // Line as P = start+t*V  where V is unit vector and t=[0..mag(end-start)]
    vector V(end - start);
    scalar magV = mag(V);
    if (magV < rootVSmall)
    {
        return 0;
    }
    V /= magV;


    // We now get the nearest intersections to start. This can either be
    // the intersection with the end plane or with the cylinder side.

    // Get the two points (expressed in t) on the end planes. This is to
    // clip any cylinder intersection against.
    scalar tPoint1;
    scalar tPoint2;

    // Maintain the two intersections with the endcaps
    scalar tNearI = vGreat;
    scalar tNearO = vGreat;
    scalar tFarI = vGreat;
    scalar tFarO = vGreat;

    // Only one cap check is required since hitting the inner radius cap is not
    // a hit
    {
        scalar s = (V & unitDir_);
        if (mag(s) > vSmall)
        {
            tPoint1 = -s1/s;
            tPoint2 = -(point2Start & unitDir_)/s;
            if (tPoint2 < tPoint1)
            {
                Swap(tPoint1, tPoint2);
            }
            // First intersection with cap is "above" or last
            // intersection is "below"
            if (tPoint1 > magV || tPoint2 < 0)
            {
                return 0;
            }

            if (tPoint1 >= 0 && tPoint1 <= magV)
            {
                const scalar rSqr = radius2(start + tPoint1*V);
                if (rSqr <= sqr(outerRadius_) && rSqr > sqr(innerRadius_))
                {
                    tNearO = tPoint1;
                }
            }
            if (tPoint2 >= 0 && tPoint2 <= magV)
            {
                const scalar rSqr = radius2(start + tPoint1*V);
                if (rSqr <= sqr(outerRadius_) && rSqr > sqr(innerRadius_))
                {
                    // Check if already have a near hit from point1
                    if (tNearO <= magV)
                    {
                        tFarO = tPoint2;
                    }
                    else
                    {
                        tNearO = tPoint2;
                    }
                }
            }
        }
        else
        {
            // Vector perpendicular to cylinder. Check for outside already done
            // above so just set tpointO to allow all.
            tPoint1 = -vGreat;
            tPoint2 = vGreat;
        }
    }

    const vector x = point1Start ^ unitDir_;
    const vector y = V ^ unitDir_;
    const scalar dI = sqr(innerRadius_);
    const scalar dO = sqr(outerRadius_);

    // Second order equation of the form a*t^2 + b*t + c
    const scalar a = (y & y);
    const scalar b = 2*(x & y);
    const scalar cI = (x & x) - dI;
    const scalar cO = (x & x) - dO;

    const scalar discI = b*b - 4*a*cI;
    const scalar discO = b*b - 4*a*cO;

    // Fully outside
    if (discO < 0)
    {
        return 0;
    }

    // Aligned with axis. Check if outsid outer radius or inside inner radius
    else if (mag(a) < rootVSmall && (cO > 0 || cI < 0))
    {
        return 0;
    }

    else if (mag(a) > rootVSmall)
    {
        // Check inner intersections, no cap intersections so fewer required checks
        if (discI < 0)
        {}
        else if (discI < rootVSmall)
        {
            // Single solution
            const scalar t = -b/(2*a);

            // With in the bounds of the line
            if (t >= 0 && t <= magV)
            {
                // Only intersection point added
                tNearI = t;
            }
        }
        else if (mag(a) > rootVSmall)
        {
            const scalar sqrtDisc = sqrt(discI);

            scalar t1 = (-b - sqrtDisc)/(2*a);
            scalar t2 = (-b + sqrtDisc)/(2*a);

            if (t2 < t1)
            {
                Swap(t1, t2);
            }

            if (t1 >= 0 && t1 <= magV)
            {
                tNearI = t1;
                if (t2 >= 0 && t2 <= magV)
                {
                    tFarI = t2;
                }
            }
            else if (t2 >= 0 && t2 <= magV)
            {
                tNearI = t2;
            }
        }

        if (discO < rootVSmall && mag(a) > rootVSmall)
        {
            // Single solution
            const scalar t = -b/(2*a);

            // Pout<< "single solution t:" << t1
            //    << " for start:" << start << " end:" << end
            //    << " c:" << c << endl;

            if (t >= 0 && t <= magV && t >= tPoint1 && t <= tPoint2)
            {
                // valid. Insert sorted.
                if (t < tNearO)
                {
                    tFarO = tNearO;
                    tNearO = t;
                }
                else if (t < tFarO)
                {
                    tFarO = t;
                }
            }
        }
        else if (mag(a) > rootVSmall)
        {
            const scalar sqrtDisc = sqrt(discO);

            scalar t1 = (-b - sqrtDisc)/(2*a);
            scalar t2 = (-b + sqrtDisc)/(2*a);
            if (t2 < t1)
            {
                Swap(t1, t2);
            }

            if (t1 >= 0 && t1 <= magV && t1 >= tPoint1 && t1 <= tPoint2)
            {
                // valid. Insert sorted.
                if (t1 < tNearO)
                {
                    tFarO = tNearO;
                    tNearO = t1;
                }
                else if (t1 < tFarO)
                {
                    tFarO = t1;
                }
            }
            if (t2 >= 0 && t2 <= magV && t2 >= tPoint1 && t2 <= tPoint2)
            {
                // valid. Insert sorted.
                if (t2 < tNearO)
                {
                    tFarO = tNearO;
                    tNearO = t2;
                }
                else if (t2 < tFarO)
                {
                    tFarO = t2;
                }
            }
            // Pout<< "two solutions t1:" << t1 << " t2:" << t2
            //    << " for start:" << start << " end:" << end
            //    << " magV:" << magV
            //    << " c:" << c << endl;
        }
    }

    label n = 0;

    // Check outside cylinder hits
    if (tNearO >= 0 && tNearO <= magV)
    {
        nearO.setPoint(start + tNearO*V);
        nearO.setHit();
        nearO.setIndex(0);
        n++;

        if (tFarO <= magV)
        {
            farO.setPoint(start + tFarO*V);
            farO.setHit();
            farO.setIndex(0);
            n++;
        }
    }
    else if (tFarO >= 0 && tFarO <= magV)
    {
        nearO.setPoint(start + tFarO*V);
        nearO.setHit();
        nearO.setIndex(0);
        n++;
    }

    // Check inner cylinder hits
    if (tNearI >= 0 && tNearI <= magV)
    {
        nearI.setPoint(start + tNearI*V);
        nearI.setHit();
        nearI.setIndex(0);
        n++;

        if (tFarI <= magV)
        {
            farI.setPoint(start + tFarI*V);
            farI.setHit();
            farI.setIndex(0);
            n++;
        }
    }
    return n;
}


Foam::boundBox Foam::searchableAnnulusCylinder::calcBounds() const
{

    // Adapted from
    // http://www.gamedev.net/community/forums
    //       /topic.asp?topic_id=338522&forum_id=20&gforum_id=0

    // Let cylinder have end points A,B and radius r,

    // Bounds in direction X (same for Y and Z) can be found as:
    // Let A.X<B.X (otherwise swap points)
    // Good approximate lowest bound is A.X-r and highest is B.X+r (precise for
    // capsule). At worst, in one direction it can be larger than needed by 2*r.

    // Accurate bounds for cylinder is
    // A.X-kx*r, B.X+kx*r
    // where
    // kx=sqrt(((A.Y-B.Y)^2+(A.Z-B.Z)^2)/((A.X-B.X)^2+(A.Y-B.Y)^2+(A.Z-B.Z)^2))

    // similar thing for Y and Z
    // (i.e.
    // ky=sqrt(((A.X-B.X)^2+(A.Z-B.Z)^2)/((A.X-B.X)^2+(A.Y-B.Y)^2+(A.Z-B.Z)^2))
    // kz=sqrt(((A.X-B.X)^2+(A.Y-B.Y)^2)/((A.X-B.X)^2+(A.Y-B.Y)^2+(A.Z-B.Z)^2))
    // )

    // How derived: geometric reasoning. Bounds of cylinder is same as for 2
    // circles centered on A and B. This sqrt thingy gives sine of angle between
    // axis and direction, used to find projection of radius.

    vector kr
    (
        sqrt(sqr(unitDir_.y()) + sqr(unitDir_.z())),
        sqrt(sqr(unitDir_.x()) + sqr(unitDir_.z())),
        sqrt(sqr(unitDir_.x()) + sqr(unitDir_.y()))
    );

    kr *= outerRadius_;

    point min = point1_ - kr;
    point max = point1_ + kr;

    min = ::Foam::min(min, point2_ - kr);
    max = ::Foam::max(max, point2_ + kr);

    return boundBox(min, max);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableAnnulusCylinder::searchableAnnulusCylinder
(
    const IOobject& io,
    const point& point1,
    const point& point2,
    const scalar innerRadius,
    const scalar outerRadius
)
:
    searchableSurface(io),
    point1_(point1),
    point2_(point2),
    magDir_(mag(point2_-point1_)),
    unitDir_((point2_-point1_)/magDir_),
    innerRadius_(innerRadius),
    outerRadius_(outerRadius)
{
    bounds() = calcBounds();
}


Foam::searchableAnnulusCylinder::searchableAnnulusCylinder
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableSurface(io),
    point1_(dict.lookup("point1")),
    point2_(dict.lookup("point2")),
    magDir_(mag(point2_-point1_)),
    unitDir_((point2_-point1_)/magDir_),
    innerRadius_(dict.lookup<scalar>("innerRadius")),
    outerRadius_(dict.lookup<scalar>("outerRadius"))
{
    bounds() = calcBounds();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchableAnnulusCylinder::~searchableAnnulusCylinder()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::searchableAnnulusCylinder::regions() const
{
    if (regions_.empty())
    {
        regions_.setSize(1);
        regions_[0] = "region0";
    }
    return regions_;
}


void Foam::searchableAnnulusCylinder::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    info.setSize(samples.size());
    forAll(samples, i)
    {
        info[i] = findNearest(samples[i], nearestDistSqr[i]);
    }
}


void Foam::searchableAnnulusCylinder::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.setSize(start.size());

    pointIndexHit b, c, d;
    forAll(start, i)
    {
        // Pick nearest intersection. If none intersected last one.
        findLineAll(start[i], end[i], info[i], b, c, d);
        if (info[i].hit())
        {}
        else if (b.hit())
        {
            info[i] = b;
        }
        else if (c.hit())
        {
            info[i] = c;
        }
        else
        {
            info[i] = d;
        }
    }
}


void Foam::searchableAnnulusCylinder::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.setSize(start.size());

    pointIndexHit b, c, d;
    forAll(start, i)
    {
        // Pick nearest intersection. If none intersected last one.
        findLineAll(start[i], end[i], info[i], b, c, d);
        if (info[i].hit())
        {}
        else if (b.hit())
        {
            info[i] = b;
        }
        else if (c.hit())
        {
            info[i] = c;
        }
        else
        {
            info[i] = d;
        }
    }
}


void Foam::searchableAnnulusCylinder::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit>>& info
) const
{
    info.setSize(start.size());

    pointIndexHit nearO, nearI, farI, farO;
    forAll(start, i)
    {
        const label n = findLineAll(start[i], end[i], nearO, nearI, farI, farO);

        info[i].setSize(n);
        label hiti = 0;
        if (nearO.hit())
        {
            info[i][hiti++] = nearO;
        }
        if (nearI.hit())
        {
            info[i][hiti++] = nearI;
        }
        if (farI.hit())
        {
            info[i][hiti++] = farI;
        }
        if (farO.hit())
        {
            info[i][hiti++] = farO;
        }
    }
}


void Foam::searchableAnnulusCylinder::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    region.setSize(info.size());
    region = 0;
}


void Foam::searchableAnnulusCylinder::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    normal.setSize(info.size());
    normal = Zero;

    forAll(info, i)
    {
        if (info[i].hit())
        {
            vector v(info[i].hitPoint() - point1_);

            // Decompose sample-point1 into normal and parallel component
            scalar parallel = (v & unitDir_);

            // Remove the parallel component and normalise
            v -= parallel*unitDir_;
            const scalar magV = mag(v);
            const scalar dI = innerRadius_ -  magV;
            const scalar dO = magV - outerRadius_;

            if (parallel <= 0)
            {
                if (dO < mag(parallel))
                {
                    // either above endcap (magV<radius) or outside but closer
                    normal[i] = -unitDir_;
                }
                else
                {
                    normal[i] = v/magV;
                }
            }
            else if (parallel <= 0.5*magDir_)
            {
                // See if endcap closer or sidewall
                if (dI >= 0 || (-dI < parallel && dI > dO))
                {
                    normal[i] = -v/magV;
                }
                else if (dO >= 0 || -dO < parallel)
                {
                    normal[i] = v/magV;
                }
                else
                {
                    // closer to endcap
                    normal[i] = -unitDir_;
                }
            }
            else if (parallel <= magDir_)
            {
                const scalar parallel2 = magDir_ - parallel;

                // See if endcap closer or sidewall
                if (dI >= 0 || (-dI < parallel2 && dI > dO))
                {
                    normal[i] = -v/magV;
                }
                else if (dO >= 0 || -dO < parallel2)
                {
                    normal[i] = v/magV;
                }
                else
                {
                    // closer to endcap
                    normal[i] = -unitDir_;
                }
            }
            else    // beyond cylinder
            {
                if (dO < parallel - magDir_)
                {
                    // either above endcap (magV<radius) or outside but closer
                    normal[i] = -unitDir_;
                }
                else
                {
                    normal[i] = v/magV;
                }
            }
        }
    }
}


void Foam::searchableAnnulusCylinder::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    volType.setSize(points.size());
    volType = volumeType::inside;

    forAll(points, pointi)
    {
        const point& pt = points[pointi];

        vector v(pt - point1_);

        // Decompose sample-point1 into normal and parallel component
        scalar parallel = v & unitDir_;

        if (parallel < 0)
        {
            // left of point1 endcap
            volType[pointi] = volumeType::outside;
        }
        else if (parallel > magDir_)
        {
            // right of point2 endcap
            volType[pointi] = volumeType::outside;
        }
        else
        {
            // Remove the parallel component
            v -= parallel*unitDir_;
            const scalar magV = mag(v);

            if (magV > outerRadius_ || magV < innerRadius_)
            {
                volType[pointi] = volumeType::outside;
            }
            else
            {
                volType[pointi] = volumeType::inside;
            }
        }
    }
}


// ************************************************************************* //
