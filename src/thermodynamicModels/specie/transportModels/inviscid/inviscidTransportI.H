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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::inviscidTransport<Thermo>::inviscidTransport
(
    const Thermo& t
)
:
    Thermo(t)
{}


template<class Thermo>
inline Foam::inviscidTransport<Thermo>::inviscidTransport
(
    const word& name,
    const inviscidTransport& ct
)
:
    Thermo(name, ct)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::inviscidTransport<Thermo>>
Foam::inviscidTransport<Thermo>::clone() const
{
    return autoPtr<inviscidTransport<Thermo>>
    (
        new inviscidTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::inviscidTransport<Thermo>>
Foam::inviscidTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<inviscidTransport<Thermo>>
    (
        new inviscidTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::inviscidTransport<Thermo>::mu
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Thermo>
inline Foam::scalar Foam::inviscidTransport<Thermo>::kappa
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::inviscidTransport<Thermo>::operator+=
(
    const inviscidTransport<Thermo>& st
)
{
    Thermo::operator+=(st);
}


template<class Thermo>
inline void Foam::inviscidTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::inviscidTransport<Thermo> Foam::operator+
(
    const inviscidTransport<Thermo>& ct1,
    const inviscidTransport<Thermo>& ct2
)
{
    return inviscidTransport<Thermo>
    (
        Thermo
        (
            static_cast<const Thermo&>(ct1)
          + static_cast<const Thermo&>(ct2)
        )
    );
}


template<class Thermo>
inline Foam::inviscidTransport<Thermo> Foam::operator*
(
    const scalar s,
    const inviscidTransport<Thermo>& ct
)
{
    return inviscidTransport<Thermo>(s*static_cast<const Thermo&>(ct));
}


// ************************************************************************* //
