/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a density based thermodynamic
                            class
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::perfectGas<Specie>::perfectGas
(
    const Specie& sp
)
:
    Specie(sp)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::perfectGas<Specie>::perfectGas
(
    const word& name,
    const perfectGas<Specie>& pf
)
:
    Specie(name, pf)
{}


template<class Specie>
inline Foam::autoPtr<Foam::perfectGas<Specie>>
Foam::perfectGas<Specie>::clone() const
{
    return autoPtr<perfectGas<Specie>>
    (
        new perfectGas<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::perfectGas<Specie>>
Foam::perfectGas<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<perfectGas<Specie>>
    (
        new perfectGas<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
Foam::scalar Foam::perfectGas<Specie>::p
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const bool limit
) const
{
    return limit ? max(pRhoT(rho, T), 0.0) : pRhoT(rho, T);
}


template<class Specie>
Foam::scalar Foam::perfectGas<Specie>::pRhoT
(
    const scalar rho,
    const scalar T
) const
{
    return rho*this->R()*T;
}


template<class Specie>
Foam::scalar Foam::perfectGas<Specie>::Gamma
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar cv
) const
{
    return this->R()/cv + 1.0;
}


template<class Specie>
Foam::scalar Foam::perfectGas<Specie>::cSqr
(
    const scalar p,
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar cv
) const
{
    return Gamma(rho, e, T, cv)*T*this->R();
}


template<class Specie>
Foam::scalar Foam::perfectGas<Specie>::dpdv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return -this->R()*T*sqr(rho);
}


template<class Specie>
Foam::scalar Foam::perfectGas<Specie>::dpde
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    NotImplemented;
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::perfectGas<Specie>::dpdT
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return rho*this->R();
}


template<class Specie>
Foam::scalar Foam::perfectGas<Specie>::E
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::perfectGas<Specie>::Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::perfectGas<Specie>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::perfectGas<Specie>::H
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::perfectGas<Specie>::CpMCv
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar CpCv
) const
{
    return this->R();
}


template<class Specie>
Foam::scalar Foam::perfectGas<Specie>::S
(
    const scalar p,
    const scalar rho,
    const scalar T
) const
{
    return 0.0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::perfectGas<Specie>::operator+=
(
    const perfectGas<Specie>& pf
)
{
    Specie::operator+=(pf);
}


template<class Specie>
inline void Foam::perfectGas<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::perfectGas<Specie> Foam::operator+
(
    const perfectGas<Specie>& pf1,
    const perfectGas<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
      + static_cast<const Specie&>(pf2)
    );

    return perfectGas<Specie>
    (
        sp
    );
}


template<class Specie>
inline Foam::perfectGas<Specie> Foam::operator*
(
    const scalar s,
    const perfectGas<Specie>& pf
)
{
    return perfectGas<Specie>
    (
        s*static_cast<const Specie&>(pf)
    );
}


template<class Specie>
inline Foam::perfectGas<Specie> Foam::operator==
(
    const perfectGas<Specie>& pf1,
    const perfectGas<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
     == static_cast<const Specie&>(pf2)
    );

    return perfectGas<Specie>
    (
        sp
    );
}


// ************************************************************************* //
