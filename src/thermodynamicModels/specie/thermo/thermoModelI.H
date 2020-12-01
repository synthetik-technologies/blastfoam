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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
inline Foam::thermoModel<ThermoType>::thermoModel
(
    const ThermoType& sp
)
:
    ThermoType(sp),
    tolerance_(1e-4),
    maxIter_(100)
{}


template<class ThermoType>
inline Foam::thermoModel<ThermoType>::thermoModel
(
    const word& name,
    const ThermoType& sp
)
:
    ThermoType(name, sp),
    tolerance_(1e-4),
    maxIter_(100)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
Foam::scalar Foam::thermoModel<ThermoType>::Gamma
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return ThermoType::Gamma(rho, e, T, ThermoType::Cv(rho, e, T));
}


template<class ThermoType>
Foam::scalar Foam::thermoModel<ThermoType>::speedOfSound
(
    const scalar& p,
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return ThermoType::speedOfSound(p, rho, e, T, ThermoType::Cv(rho, e, T));
}


template<class ThermoType>
Foam::scalar Foam::thermoModel<ThermoType>::CpByCv
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return ThermoType::Cp(rho, e, T)/ThermoType::Cv(rho, e, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class ThermoType>
inline void Foam::thermoModel<ThermoType>::operator+=
(
    const thermoModel<ThermoType>& st
)
{
    ThermoType::operator+=(st);
}


template<class ThermoType>
inline void Foam::thermoModel<ThermoType>::operator*=(const scalar s)
{
    ThermoType::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class ThermoType>
inline Foam::thermoModel<ThermoType> Foam::operator+
(
    const thermoModel<ThermoType>& st1,
    const thermoModel<ThermoType>& st2
)
{
    return thermoModel<ThermoType>
    (
        static_cast<const ThermoType&>(st1)
      + static_cast<const ThermoType&>(st2)
    );
}


template<class ThermoType>
inline Foam::thermoModel<ThermoType> Foam::operator*
(
    const scalar s,
    const thermoModel<ThermoType>& st
)
{
    return thermoModel<ThermoType>
    (
        s*static_cast<const ThermoType&>(st)
    );
}


template<class ThermoType>
inline Foam::thermoModel<ThermoType> Foam::operator==
(
    const thermoModel<ThermoType>& st1,
    const thermoModel<ThermoType>& st2
)
{
    return thermoModel<ThermoType>
    (
        static_cast<const ThermoType&>(st1)
     == static_cast<const ThermoType&>(st2)
    );
}


// ************************************************************************* //
