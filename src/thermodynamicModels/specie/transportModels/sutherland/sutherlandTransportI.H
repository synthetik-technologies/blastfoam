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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
inline void Foam::sutherlandTransport<Thermo>::calcCoeffs
(
    const scalar mu1, const scalar Smu1,
    const scalar kappa1, const scalar Sk1,
    const scalar mu2, const scalar Smu2,
    const scalar kappa2, const scalar Sk2
)
{
    scalar rootSmu1 = sqrt(Smu1);
    scalar rootSmu2 = sqrt(Smu2);
    scalar mu1rootSmu2 = mu1*rootSmu2;
    scalar mu2rootSmu1 = mu2*rootSmu1;

    Smu_ = (mu2rootSmu1 - mu1rootSmu2)/(mu1rootSmu2/Smu1 - mu2rootSmu1/Smu2);
    Amu_ = mu1*(1.0 + Smu_/Smu1)/rootSmu1;

    scalar rootSk1 = sqrt(Sk1);
    scalar rootSk2 = sqrt(Sk2);
    scalar kappa1rootSk2 = kappa1*rootSk2;
    scalar kappa2rootSk1 = kappa2*rootSk1;

    Sk_ =
        (kappa2rootSk1 - kappa1rootSk2)
       /(kappa1rootSk2/Sk1 - kappa2rootSk1/Sk2);
    Ak_ = kappa1*(1.0 + Sk_/Sk1)/rootSk1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::sutherlandTransport<Thermo>::sutherlandTransport
(
    const Thermo& t,
    const scalar Amu,
    const scalar Smu,
    const scalar Ak,
    const scalar Sk
)
:
    Thermo(t),
    Amu_(Amu),
    Smu_(Smu),
    Ak_(Ak),
    Sk_(Sk)
{}


template<class Thermo>
inline Foam::sutherlandTransport<Thermo>::sutherlandTransport
(
    const Thermo& t,
    const scalar mu1, const scalar Smu1, const scalar kappa1, const scalar Sk1,
    const scalar mu2, const scalar Smu2, const scalar kappa2, const scalar Sk2
)
:
    Thermo(t)
{
    calcCoeffs(mu1, Smu1, kappa1, Sk1, mu2, Smu2, kappa2, Sk2);
}


template<class Thermo>
inline Foam::sutherlandTransport<Thermo>::sutherlandTransport
(
    const word& name,
    const sutherlandTransport& st
)
:
    Thermo(name, st),
    Amu_(st.Amu_),
    Smu_(st.Smu_),
    Ak_(st.Ak_),
    Sk_(st.Sk_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::sutherlandTransport<Thermo>>
Foam::sutherlandTransport<Thermo>::clone() const
{
    return autoPtr<sutherlandTransport<Thermo>>
    (
        new sutherlandTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::sutherlandTransport<Thermo>>
Foam::sutherlandTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<sutherlandTransport<Thermo>>
    (
        new sutherlandTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::sutherlandTransport<Thermo>::mu
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Amu_*::sqrt(T)/(1.0 + Smu_/T);
}


template<class Thermo>
inline Foam::scalar Foam::sutherlandTransport<Thermo>::kappa
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (Ak_ < 0)
    {
        scalar Cv_ = this->Cv(rho, e, T);
        return mu(rho, e, T)*Cv_*(1.32 + 1.77*this->R()/Cv_);
    }
    else
    {
        return Ak_*::sqrt(T)/(1.0 + Sk_/T);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::sutherlandTransport<Thermo>::operator+=
(
    const sutherlandTransport<Thermo>& st
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(st);

    if ((Ak_ < 0) != (st.Ak_ < 0))
    {
        FatalErrorInFunction
            << "Thermal conductivity not specified for all species" << endl
            << exit(FatalError);
    }

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        scalar Y2 = st.Y()/this->Y();

        Amu_ = Y1*Amu_ + Y2*st.Amu_;
        Smu_ = Y1*Smu_ + Y2*st.Smu_;

        Ak_ = Y1*Ak_ + Y2*st.Ak_;
        Sk_ = Y1*Sk_ + Y2*st.Sk_;
    }
}


template<class Thermo>
inline void Foam::sutherlandTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::sutherlandTransport<Thermo> Foam::operator+
(
    const sutherlandTransport<Thermo>& st1,
    const sutherlandTransport<Thermo>& st2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );

    if ((st1.Ak_ < 0) != (st2.Ak_ < 0))
    {
        FatalErrorInFunction
            << "Thermal conductivity not specified for all species" << endl
            << exit(FatalError);
    }

    if (mag(t.Y()) < small)
    {
        return sutherlandTransport<Thermo>
        (
            t,
            0,
            st1.Amu_,
            st1.Smu_,
            st1.Ak_,
            st1.Sk_
        );
    }
    else
    {
        scalar Y1 = st1.Y()/t.Y();
        scalar Y2 = st2.Y()/t.Y();

        return sutherlandTransport<Thermo>
        (
            t,
            Y1*st1.Amu_ + Y2*st2.Amu_,
            Y1*st1.Smu_ + Y2*st2.Smu_,
            Y1*st1.Ak_ + Y2*st2.Ak_,
            Y1*st1.Sk_ + Y2*st2.Sk_
        );
    }
}


template<class Thermo>
inline Foam::sutherlandTransport<Thermo> Foam::operator*
(
    const scalar s,
    const sutherlandTransport<Thermo>& st
)
{
    return sutherlandTransport<Thermo>
    (
        s*static_cast<const Thermo&>(st),
        st.Amu_,
        st.Smu_,
        st.Ak_,
        st.Sk_
    );
}


// ************************************************************************* //
