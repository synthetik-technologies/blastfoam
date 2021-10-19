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

#include "IOstreams.H"
#include "pointHit.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Point, class PointRef>
inline triangleShape<Point, PointRef>::triangleShape
<Point, PointRef>::triangleShape
(
    const Point& a,
    const Point& b,
    const Point& c
)
:
    triangle<Point, PointRef>(a, b, c)
{}


template<class Point, class PointRef>
inline triangleShape<Point, PointRef>::triangleShape(Istream& is)
:
    triangle<Point, PointRef>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Point, class PointRef>
inline vector triangleShape<Point, PointRef>::gradNi
(
    const label i
) const
{
    vector gradN = vector::zero;

    switch(i)
    {
        case 0:
            gradN = -0.5*((this->c() - this->b())^this->normal())/sqr(this->mag());
            break;
        case 1:
            gradN = -0.5*((this->a() - this->c())^this->normal())/sqr(this->mag());
            break;
        case 2:
            gradN = -0.5*((this->b() - this->a())^this->normal())/sqr(this->mag());
            break;
        default:
            FatalErrorInFunction
                << "Index of triangle shape function is out of range."
                << abort(FatalError);
    }

    return gradN;
}


template<class Point, class PointRef>
inline scalar triangleShape<Point, PointRef>::Ni
(
    const label i,
    const Point& p
) const
{
    vector gradN = gradNi(i);
    scalar N = 0.0;

    switch(i)
    {
        case 0:
            N = (gradN & (p - this->b()));
            break;
        case 1:
            N = (gradN & (p - this->c()));
            break;
        case 2:
            N = (gradN & (p - this->a()));
            break;
        default:
            FatalErrorInFunction
                << "Index of triangle shape function is out of range."
                << abort(FatalError);
    }

    return N;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Point, class PointRef>
inline Istream& operator>>(Istream& is, triangleShape<Point, PointRef>& t)
{
    // Read beginning of triangle point pair
    is.readBegin("triangle");

    is >> t.a_ >> t.b_ >> t.c_;

    // Read end of triangle point pair
    is.readEnd("triangle");

    // Check state of Ostream
    is.check("Istream& operator>>(Istream&, triangle&)");

    return is;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

