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

inline rBlastSpecie::rBlastSpecie
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


inline rBlastSpecie::rBlastSpecie
(
    const scalar Y,
    const scalar molWeight
)
:
    Y_(Y),
    molWeight_(molWeight)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline rBlastSpecie::rBlastSpecie(const word& name, const rBlastSpecie& st)
:
    specieName_(name),
    Y_(st.Y_),
    molWeight_(st.molWeight_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const word& rBlastSpecie::specieName() const
{
    return specieName_;
}


inline scalar rBlastSpecie::W() const
{
    return molWeight_;
}


inline scalar rBlastSpecie::Y() const
{
    return Y_;
}


inline scalar rBlastSpecie::R() const
{
    return RR/molWeight_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// inline void rBlastSpecie::operator=(const rBlastSpecie& st)
// {
//     //name_ = st.name_;
//     Y_ = st.Y_;
//     molWeight_ = st.molWeight_;
// }


inline void rBlastSpecie::operator+=(const rBlastSpecie& st)
{
    const scalar sumY = Y_ + st.Y_;
    if (mag(sumY) > small)
    {
        molWeight_ = sumY/(Y_/molWeight_ + st.Y_/st.molWeight_);
    }

    Y_ = sumY;
}


inline void rBlastSpecie::operator*=(const scalar s)
{
    Y_ *= s;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline rBlastSpecie operator+(const rBlastSpecie& st1, const rBlastSpecie& st2)
{
    const scalar sumY = max(st1.Y_ + st2.Y_, small);

    if (mag(sumY) > small)
    {
        return rBlastSpecie
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


inline rBlastSpecie operator*(const scalar s, const rBlastSpecie& st)
{
    return rBlastSpecie
    (
        s*st.Y_,
        st.molWeight_
    );
}


inline rBlastSpecie operator==(const rBlastSpecie& st1, const rBlastSpecie& st2)
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

    return rBlastSpecie(diffY, molWeight);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
