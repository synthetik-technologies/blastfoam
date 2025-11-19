/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "specieBlast.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::keyesTransport<Thermo>::keyesTransport
(
    const Thermo& t,
    const scalar mu0,
    const scalar B,
    const scalar C,
    const scalar rPr
)
:
    Thermo(t),
    mu0_(mu0),
    B_(B),
    C_(C),
    rPr_(rPr)
{}


template<class Thermo>
inline Foam::keyesTransport<Thermo>::keyesTransport
(
    const word& name,
    const keyesTransport& st
)
:
    Thermo(name, st),
    mu0_(st.mu0_),
    B_(st.B_),
    C_(st.C_),
    rPr_(st.rPr_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::keyesTransport<Thermo>>
Foam::keyesTransport<Thermo>::clone() const
{
    return autoPtr<keyesTransport<Thermo>>
    (
        new keyesTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::keyesTransport<Thermo>>
Foam::keyesTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<keyesTransport<Thermo>>
    (
        new keyesTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::keyesTransport<Thermo>::mu
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return mu0_*::sqrt(T)/(1.0 + (B_/T)*pow(10, -C_/T));
}


template<class Thermo>
inline Foam::scalar Foam::keyesTransport<Thermo>::kappa
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    // if (rPr_ > 0)
    {
        return mu(rho, e, T)*this->Cp(rho, e, T)*rPr_;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::keyesTransport<Thermo>::operator+=
(
    const keyesTransport<Thermo>& st
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(st);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        scalar Y2 = st.Y()/this->Y();

        mu0_ = Y1*mu0_ + Y2*st.mu0_;
        B_ = Y1*B_ + Y2*st.B_;
        C_ = Y1*C_ + Y2*st.C_;
        rPr_ = 1.0/(Y1/rPr_ + Y2/st.rPr_);
    }
}


template<class Thermo>
inline void Foam::keyesTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::keyesTransport<Thermo> Foam::operator+
(
    const keyesTransport<Thermo>& st1,
    const keyesTransport<Thermo>& st2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );

    if (mag(t.Y()) < small)
    {
        return keyesTransport<Thermo>
        (
            t,
            0,
            st1.mu0_,
            st1.B_,
            st1.C_,
            st1.rPr_
        );
    }
    else
    {
        scalar Y1 = st1.Y()/t.Y();
        scalar Y2 = st2.Y()/t.Y();

        return keyesTransport<Thermo>
        (
            t,
            Y1*st1.mu0_ + Y2*st2.mu0_,
            Y1*st1.B_ + Y2*st2.B_,
            Y1*st1.C_ + Y2*st2.C_,
            1.0/(Y1/st1.rPr_ + Y2/st2.rPr_)
        );
    }
}


template<class Thermo>
inline Foam::keyesTransport<Thermo> Foam::operator*
(
    const scalar s,
    const keyesTransport<Thermo>& st
)
{
    return keyesTransport<Thermo>
    (
        s*static_cast<const Thermo&>(st),
        st.mu0_,
        st.B_,
        st.C_,
        st.rPr_
    );
}


// ************************************************************************* //
