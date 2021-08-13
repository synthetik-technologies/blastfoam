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
inline Foam::constAnIsoSolidTransport<Thermo>::constAnIsoSolidTransport
(
    const Thermo& t,
    const vector kappa
)
:
    Thermo(t),
    kappa_(kappa)
{}


template<class Thermo>
inline Foam::constAnIsoSolidTransport<Thermo>::constAnIsoSolidTransport
(
    const word& name,
    const constAnIsoSolidTransport& ct
)
:
    Thermo(name, ct),
    kappa_(ct.kappa_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::constAnIsoSolidTransport<Thermo>::kappa
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return mag(kappa_);
}

template<class Thermo>
inline Foam::vector Foam::constAnIsoSolidTransport<Thermo>::Kappa
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return kappa_;
}


template<class Thermo>
inline Foam::scalar Foam::constAnIsoSolidTransport<Thermo>::mu
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
inline Foam::scalar Foam::constAnIsoSolidTransport<Thermo>::alphah
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return mag(kappa_)/this->Cp(rho, e, T);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::constAnIsoSolidTransport<Thermo>::operator+=
(
    const constAnIsoSolidTransport<Thermo>& ct
)
{
    scalar Y1 = this->Y();

    Y1 /= this->Y();
    scalar Y2 = ct.Y()/this->Y();

    kappa_ = Y1*kappa_ + Y2*ct.kappa_;
}


template<class Thermo>
inline void Foam::constAnIsoSolidTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::constAnIsoSolidTransport<Thermo> Foam::operator+
(
    const constAnIsoSolidTransport<Thermo>& ct1,
    const constAnIsoSolidTransport<Thermo>& ct2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ct1) + static_cast<const Thermo&>(ct2)
    );

    if (mag(t.Y()) < small)
    {
        return constAnIsoSolidTransport<Thermo>
        (
            t,
            ct1.kappa_
        );
    }
    else
    {
        scalar Y1 = ct1.Y()/t.Y();
        scalar Y2 = ct2.Y()/t.Y();

        return constAnIsoSolidTransport<Thermo>
        (
            t,
            Y1*ct1.kappa_ + Y2*ct2.kappa_
        );
    }
}


template<class Thermo>
inline Foam::constAnIsoSolidTransport<Thermo> Foam::operator*
(
    const scalar s,
    const constAnIsoSolidTransport<Thermo>& ct
)
{
    return constAnIsoSolidTransport<Thermo>
    (
        s*static_cast<const Thermo&>(ct),
        ct.kappa_
    );
}

// ************************************************************************* //
