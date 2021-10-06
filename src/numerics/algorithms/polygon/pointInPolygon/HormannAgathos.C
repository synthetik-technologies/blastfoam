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

Class
  HormannAgathos

Description
    Implements the point in polygon problem using the winding number technique
    presented in the paper:
        "The point in polygon problem for arbitrary polygons",
        Kai Hormann, Alexander Agathos, 2001

Author
    Martin Beaudoin, Hydro-Quebec, (2008)

\*----------------------------------------------------------------------------*/

#include "HormannAgathos.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HormannAgathos::HormannAgathos
(
    const List<point2D>& P,
    const scalar& distTol
)
:
    P_(P),
    distTol_(distTol)
{
    evaluateEpsilon();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// "The point in polygon problem for arbitrary polygons",
// Kai Hormann , Alexander Agathos
// We implement the algorthm #7 of the paper, the "Efficient boundary algorithm"
Foam::HormannAgathos::inOutClassification
Foam::HormannAgathos::evaluate(const point2D& R) const
{
    inOutClassification retVal = POINT_OUTSIDE;

    label w = 0;
    scalar determinantTol = 8.0*sqr(epsilon_);

    // For this part of the algorithm, we first try to determine
    // if we are on a vertex or on a edge.
    // We need to make these comparisons using a tolerance because we are
    // manipulating floating point values.
    // The best would be to redefine point2D with a redefinition of operator =
    // and operator > to take care of comparison with tolerance...
    // I chose instead to simply implement quick inline member functions
    // to do the job.

    label indexI   = P_.size()-1;
    label indexIp1 = indexI;

    // if P[0].y == R.y) &&  P[0].x == R.x
    if
    (
        equalWithTol(P_[indexI].y(), R.y())
     && equalWithTol(P_[indexI].x(), R.x())
    )
    {
        retVal = POINT_ON_VERTEX;
    }
    else
    {
        // for i=0 to n - 1
        for (label i=0; i < P_.size(); i++)
        {
            indexI = indexIp1;
            indexIp1 = i;

            // if P[i+1].y == R.y
            if(equalWithTol(P_[indexIp1].y(), R.y()))
            {
                // if P[i+1].x == R.x
                if(equalWithTol(P_[indexIp1].x(), R.x()))
                {
                    retVal = POINT_ON_VERTEX;
                    break;
                }
                // else if P[i].y == R.y && (P[i+i].x > R.x) == (P[i].x < R.x)
                else if
                (
                    equalWithTol(P_[indexI].y(), R.y())
                 && (
                        greaterWithTol(P_[indexIp1].x(), R.x())
                     == smallerWithTol(P_[indexI].x(), R.x())
                    )
                )
                {
                    retVal = POINT_ON_EDGE;
                    break;
                }

            }

            // From here, I am not sure if I still need to use the
            // "tolerance aware" operators ==, >, <, et....
            // Need to run some more validation later on over specific
            // cases where we go play near the edges...
            // I have left the "non tolerance" aware instructions in
            // comments, jut in case.
            // MB 10/03/2008
            // if crossing
            //if ((P_[indexI].y() < R.y()) != (P_[indexIp1].y() < R.y()) )
            if
            (
                smallerWithTol(P_[indexI].y(), R.y())
             != smallerWithTol(P_[indexIp1].y(), R.y())
            )
            {
                //if(P_[indexI].x() >= R.x())
                if (greaterOrEqualWithTol(P_[indexI].x(), R.x()))
                {
                    //if(P_[indexIp1].x() > R.x())
                    if (greaterWithTol(P_[indexIp1].x(), R.x()))
                    {
                        // Modify w
                        if (greaterWithTol(P_[indexIp1].y(), P_[indexI].y()))
                        {
                            w++;
                        }
                        else
                        {
                            w--;
                        }

                    }
                    else
                    {
                        // Compute det(i)
                        scalar detI =
                            (P_[indexI].x() - R.x()) * (P_[indexIp1].y() - R.y())
                          - (P_[indexIp1].x() - R.x()) * (P_[indexI].y() - R.y());

                        // VSMALL is too strict... we need to take
                        // into account that we tolerate an epsilon
                        // error all along so in the worst case, if
                        // the point R is separated by only epsilon
                        // from an edge, the smallest determinant
                        // should be of the order of
                        // 2*(2*epsilon*2epsilon) = 8*epsilon^2
                        if(mag(detI) < determinantTol)
                        {
                            retVal = POINT_ON_EDGE;
                            break;
                        }

                        // if right_crossing
                        //if ( (detI > 0.0) == (P_[indexIp1].y() > P_[indexI].y()) )
                        if
                        (
                            (detI > 0.0)
                         == (greaterWithTol(P_[indexIp1].y(), P_[indexI].y()))
                        )
                        {
                            // Modify w
                            if (greaterWithTol(P_[indexIp1].y(), P_[indexI].y()))
                            {
                                w++;
                            }
                            else
                            {
                                w--;
                            }
                        }
                    }
                }
                //else if (P_[indexIp1].x() > R.x())
                else if (greaterWithTol(P_[indexIp1].x(), R.x()))
                {
                    // Compute det(i)
                    scalar detI =
                        (P_[indexI].x() - R.x()) * (P_[indexIp1].y() - R.y())
                      - (P_[indexIp1].x() - R.x()) * (P_[indexI].y() - R.y());

                    // VSMALL is too strict... we need to take into
                    // account that we tolerate an epsilon error all
                    // along // so in the worst case, if the point R
                    // is separated by only epsilon from an edge, //
                    // the smallest determinant should be of the order
                    // of 2*(2*epsilon*2epsilon) = 8*epsilon^2
                    if (mag(detI) < determinantTol)
                    {
                        retVal = POINT_ON_EDGE;
                        break;
                    }

                    // if right_crossing
                    //if( (detI > 0.0) == (P_[indexIp1].y() > P_[indexI].y()) )
                    if
                    (
                        (detI > 0.0)
                     == (greaterWithTol(P_[indexIp1].y(), P_[indexI].y()))
                    )
                    {
                        // Modify w
                        if (greaterWithTol(P_[indexIp1].y(), P_[indexI].y()))
                        {
                            w++;
                        }
                        else
                        {
                            w--;
                        }
                    }
                }
            }
        }
    }

    if (retVal == POINT_OUTSIDE && w != 0)
    {
        retVal = POINT_INSIDE;
    }

    return retVal;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::HormannAgathos::evaluateEpsilon()
{
    // Compute the 2D epsilon distance to detect if a point is close
    // to another We take the length of the shortest edge of the
    // polygon P_, and we multiply by distTol_

    scalar minDist2 = GREAT;
    label indexIp1 = P_.size() - 1;
    label indexI;

    scalar dist2;
    // for i=0 to n - 1
    for (label i = 0; i < P_.size(); i++)
    {
        indexI = indexIp1;
        indexIp1 = i;

        dist2 = sqr(P_[indexIp1].y() - P_[indexI].y())
            + sqr(P_[indexIp1].x() - P_[indexI].x());

        if (dist2 < minDist2)
        {
            minDist2 = dist2;
        }
    }

    epsilon_ =  distTol_ * sqrt(minDist2);
}


// ************************************************************************* //
