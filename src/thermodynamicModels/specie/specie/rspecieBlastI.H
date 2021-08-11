/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Rename class
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline rspecieBlast::rspecieBlast
(
    const word& name,
    const scalar Y,
    const scalar molWeight
)
:
    specieName_(name),
    Y_(Y),
    molWeight_(molWeight)
{}


inline rspecieBlast::rspecieBlast
(
    const scalar Y,
    const scalar molWeight
)
:
    Y_(Y),
    molWeight_(molWeight)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline rspecieBlast::rspecieBlast(const word& name, const rspecieBlast& st)
:
    specieName_(name),
    Y_(st.Y_),
    molWeight_(st.molWeight_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const word& rspecieBlast::specieName() const
{
    return specieName_;
}


inline scalar rspecieBlast::W() const
{
    return molWeight_;
}


inline scalar rspecieBlast::Y() const
{
    return Y_;
}


inline scalar rspecieBlast::R() const
{
    return RR/molWeight_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// inline void rspecieBlast::operator=(const rspecieBlast& st)
// {
//     //name_ = st.name_;
//     Y_ = st.Y_;
//     molWeight_ = st.molWeight_;
// }


inline void rspecieBlast::operator+=(const rspecieBlast& st)
{
    const scalar sumY = Y_ + st.Y_;
    if (mag(sumY) > small)
    {
        molWeight_ = sumY/(Y_/molWeight_ + st.Y_/st.molWeight_);
    }

    Y_ = sumY;
}


inline void rspecieBlast::operator*=(const scalar s)
{
    Y_ *= s;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline rspecieBlast operator+(const rspecieBlast& st1, const rspecieBlast& st2)
{
    const scalar sumY = max(st1.Y_ + st2.Y_, small);

    if (mag(sumY) > small)
    {
        return rspecieBlast
        (
            sumY,
            sumY/(st1.Y_/st1.molWeight_ + st2.Y_/st2.molWeight_)
        );
    }
    else
    {
        return st1;
    }
}


inline rspecieBlast operator*(const scalar s, const rspecieBlast& st)
{
    return rspecieBlast
    (
        s*st.Y_,
        st.molWeight_
    );
}


inline rspecieBlast operator==(const rspecieBlast& st1, const rspecieBlast& st2)
{
    scalar diffY = st2.Y_ - st1.Y_;
    if (mag(diffY) < small)
    {
        diffY = small;
    }

    const scalar diffRW =
        st2.Y_/st2.molWeight_ - st1.Y_/st1.molWeight_;

    scalar molWeight = great;
    if (mag(diffRW) > small)
    {
        molWeight = diffY/diffRW;
    }

    return rspecieBlast(diffY, molWeight);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
