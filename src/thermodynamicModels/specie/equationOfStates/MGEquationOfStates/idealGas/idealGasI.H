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
inline Foam::idealGas<Specie>::idealGas
(
    const Specie& sp,
    const scalar gamma
)
:
    Specie(sp),
    gamma_(gamma)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::idealGas<Specie>::idealGas
(
    const word& name,
    const idealGas<Specie>& pf
)
:
    Specie(name, pf),
    gamma_(pf.gamma_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::idealGas<Specie>>
Foam::idealGas<Specie>::clone() const
{
    return autoPtr<idealGas<Specie>>
    (
        new idealGas<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::idealGas<Specie>>
Foam::idealGas<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<idealGas<Specie>>
    (
        new idealGas<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
Foam::scalar Foam::idealGas<Specie>::p
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const bool limit
) const
{
    return
        limit
      ? max((gamma_ - 1.0)*rho*e, 0.0)
      : (gamma_ - 1.0)*rho*e;
}


template<class Specie>
Foam::scalar Foam::idealGas<Specie>::Gamma
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar Cv
) const
{
    return gamma_;
}


template<class Specie>
Foam::scalar Foam::idealGas<Specie>::cSqr
(
    const scalar p,
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar Cv
) const
{
    return p*gamma_/max(rho, 1e-10);
}


template<class Specie>
Foam::scalar Foam::idealGas<Specie>::dpdv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return -sqr(rho)*(gamma_ - 1.0)*e;
}


template<class Specie>
Foam::scalar Foam::idealGas<Specie>::dpde
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return (gamma_ - 1.0)*rho;
}


template<class Specie>
Foam::scalar Foam::idealGas<Specie>::dpdT
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
Foam::scalar Foam::idealGas<Specie>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::idealGas<Specie>::Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::idealGas<Specie>::CpMCv
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar CpCv
) const
{
    return this->eBased_ ? (gamma_ - 1.0)*CpCv : CpCv*(gamma_ - 1.0)/gamma_;
}


template<class Specie>
Foam::scalar Foam::idealGas<Specie>::E
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::idealGas<Specie>::H
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::idealGas<Specie>::S
(
    const scalar p,
    const scalar rho,
    const scalar T
) const
{
    return -this->R()*log(p/Foam::constant::thermodynamic::Pstd);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::idealGas<Specie>::operator+=
(
    const idealGas<Specie>& pf
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(pf);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = pf.Y()/this->Y();

        gamma_ = Y1*gamma_ + Y2*pf.gamma_;
    }
}


template<class Specie>
inline void Foam::idealGas<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::idealGas<Specie> Foam::operator+
(
    const idealGas<Specie>& pf1,
    const idealGas<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
      + static_cast<const Specie&>(pf2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return idealGas<Specie>
        (
            sp,
            pf1.gamma_
        );
    }
    else
    {
        const scalar Y1 = pf1.Y()/sp.Y();
        const scalar Y2 = pf2.Y()/sp.Y();

        return idealGas<Specie>
        (
            sp,
            Y1*pf1.gamma_ + Y2*pf2.gamma_
        );
    }
}


template<class Specie>
inline Foam::idealGas<Specie> Foam::operator*
(
    const scalar s,
    const idealGas<Specie>& pf
)
{
    return idealGas<Specie>
    (
        s*static_cast<const Specie&>(pf),
        pf.gamma_
    );
}


template<class Specie>
inline Foam::idealGas<Specie> Foam::operator==
(
    const idealGas<Specie>& pf1,
    const idealGas<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
     == static_cast<const Specie&>(pf2)
    );

    const scalar Y1 = pf1.Y()/sp.Y();
    const scalar Y2 = pf2.Y()/sp.Y();

    return idealGas<Specie>
    (
        sp,
        Y2*pf2.gamma_ - Y1*pf1.gamma_
    );
}


// ************************************************************************* //
