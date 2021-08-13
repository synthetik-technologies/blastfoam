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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::constIsoSolidTransport<Thermo>::constIsoSolidTransport
(
    const Thermo& t,
    const scalar kappa
)
:
    Thermo(t),
    kappa_(kappa)
{}


template<class Thermo>
inline Foam::constIsoSolidTransport<Thermo>::constIsoSolidTransport
(
    const word& name,
    const constIsoSolidTransport& ct
)
:
    Thermo(name, ct),
    kappa_(ct.kappa_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::constIsoSolidTransport<Thermo>::kappa
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return kappa_;
}

template<class Thermo>
inline Foam::vector Foam::constIsoSolidTransport<Thermo>::Kappa
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return vector(kappa_, kappa_, kappa_);
}


template<class Thermo>
inline Foam::scalar Foam::constIsoSolidTransport<Thermo>::mu
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    NotImplemented;
    return scalar(0);
}


template<class Thermo>
inline Foam::scalar Foam::constIsoSolidTransport<Thermo>::alphah
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return kappa_/this->Cp(rho, e, T);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::constIsoSolidTransport<Thermo>::operator+=
(
    const constIsoSolidTransport<Thermo>& ct
)
{
    scalar Y1 = this->Y();

    Y1 /= this->Y();
    scalar Y2 = ct.Y()/this->Y();

    kappa_ = Y1*kappa_ + Y2*ct.kappa_;
}


template<class Thermo>
inline void Foam::constIsoSolidTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::constIsoSolidTransport<Thermo> Foam::operator+
(
    const constIsoSolidTransport<Thermo>& ct1,
    const constIsoSolidTransport<Thermo>& ct2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ct1) + static_cast<const Thermo&>(ct2)
    );

    if (mag(t.Y()) < small)
    {
        return constIsoSolidTransport<Thermo>
        (
            t,
            ct1.kappa_
        );
    }
    else
    {
        scalar Y1 = ct1.Y()/t.Y();
        scalar Y2 = ct2.Y()/t.Y();

        return constIsoSolidTransport<Thermo>
        (
            t,
            Y1*ct1.kappa_ + Y2*ct2.kappa_
        );
    }
}


template<class Thermo>
inline Foam::constIsoSolidTransport<Thermo> Foam::operator*
(
    const scalar s,
    const constIsoSolidTransport<Thermo>& ct
)
{
    return constIsoSolidTransport<Thermo>
    (
        s*static_cast<const Thermo&>(ct),
        ct.kappa_
    );
}

// ************************************************************************* //
