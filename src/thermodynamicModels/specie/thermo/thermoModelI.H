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


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::G
(
    const scalar p,
    const scalar rho,
    const scalar T
) const
{
    scalar e = this->Es(rho, 1000.0, T);
    return this->Ha(rho, e, T) - T*this->S(p, rho, e, T);
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::A
(
    const scalar p,
    const scalar rho,
    const scalar T
) const
{
    scalar e = this->Es(rho, 1000.0, T);
    return this->Ea(rho, e, T) - T*this->S(p, rho, e, T);
}


template<class ThermoType>
inline Foam::scalar
Foam::thermoModel<ThermoType>::K
(
    const scalar p,
    const scalar rho,
    const scalar T
) const
{
    scalar arg = -this->Y()*this->G(Pstd, rho, T)/(RR*T);

    if (arg < 600)
    {
        return exp(arg);
    }
    else
    {
        return rootVGreat;
    }
}


template<class ThermoType>
inline Foam::scalar
Foam::thermoModel<ThermoType>::Kp
(
    const scalar p,
    const scalar rho,
    const scalar T
) const
{
    return K(p, rho, T);
}


template<class ThermoType>
inline Foam::scalar
Foam::thermoModel<ThermoType>::Kc
(
    const scalar p,
    const scalar rho,
    const scalar T
) const
{
    const scalar nm = this->Y()/this->W();

    if (equal(nm, small))
    {
        return Kp(p, rho, T);
    }
    else
    {
        return Kp(p, rho, T)*pow(Pstd/(RR*T), nm);
    }
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::Kx
(
    const scalar p,
    const scalar rho,
    const scalar T
) const
{
    const scalar nm = this->Y()/this->W();

    if (equal(nm, small))
    {
        return Kp(p, rho, T);
    }
    else
    {
        return Kp(p, rho, T)*pow(Pstd/p, nm);
    }
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::Kn
(
    const scalar p,
    const scalar rho,
    const scalar T,
    const scalar n
) const
{
    const scalar nm = this->Y()/this->W();

    if (equal(nm, small))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, rho, T)*pow(n*Pstd/p, nm);
    }
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return this->Cp(rho, e, T)*this->W();
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::ha
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return this->Ha(rho, e, T)*this->W();
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::hs
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return this->Hs(rho, e, T)*this->W();
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::hc() const
{
    return this->Hf()*this->W();
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::s
(
    const scalar p,
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return this->S(p, rho, e, T)*this->W();
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::he
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return this->es(rho, e, T);
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return this->Cv(rho, e, T)*this->W();
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::es
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return this->Es(rho, e, T)*this->W();
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::ea
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return this->Ea(rho, e, T)*this->W();
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::g
(
    const scalar p,
    const scalar rho,
    const scalar T
) const
{
    return this->G(p, rho, T)*this->W();
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::a
(
    const scalar p,
    const scalar rho,
    const scalar T
) const
{
    return this->A(p, rho, T)*this->W();
}



template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::dKcdTbyKc
(
    const scalar p,
    const scalar T
) const
{
    const scalar nm = this->Y()/this->W();

    if (equal(nm, small))
    {
        return -this->dGdT(Pstd, T)*this->Y()/RR;
    }
    else
    {
        return -(nm/T + this->dGdT(Pstd, T)*this->Y()/RR);
    }
}

template<class ThermoType>
inline Foam::scalar
Foam::thermoModel<ThermoType>::dcpdT(const scalar p, const scalar T) const
{
    return this->dCpdT(p, T)*this->W();
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
