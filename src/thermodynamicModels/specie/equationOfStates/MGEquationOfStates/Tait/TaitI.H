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

#include "constants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::Tait<Specie>::Tait
(
    const Specie& sp,
    const scalar a,
    const scalar b,
    const scalar gamma
)
:
    Specie(sp),
    a_(a),
    b_(b),
    gamma_(gamma)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::Tait<Specie>::Tait
(
    const word& name,
    const Tait<Specie>& pf
)
:
    Specie(name, pf),
    a_(pf.a_),
    b_(pf.b_),
    gamma_(pf.gamma_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::Tait<Specie>>
Foam::Tait<Specie>::clone() const
{
    return autoPtr<Tait<Specie>>
    (
        new Tait<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::Tait<Specie>>
Foam::Tait<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<Tait<Specie>>
    (
        new Tait<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
Foam::scalar Foam::Tait<Specie>::Pi
(
    const scalar& rho,
    const scalar& e
) const
{
    return gamma_*(b_ - a_);
}


template<class Specie>
Foam::scalar Foam::Tait<Specie>::delta
(
    const scalar& p,
    const scalar& rho,
    const scalar& e
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::Tait<Specie>::Gamma
(
    const scalar& rho,
    const scalar& e
) const
{
    return gamma_;
}


template<class Specie>
Foam::scalar Foam::Tait<Specie>::dpdv
(
    const scalar& rho,
    const scalar& e
) const
{
    return -sqr(rho)*(gamma_ - 1.0)*e;
}


template<class Specie>
Foam::scalar Foam::Tait<Specie>::dpde
(
    const scalar& rho,
    const scalar& e
) const
{
    return (gamma_ - 1.0)*rho;
}


template<class Specie>
Foam::scalar Foam::Tait<Specie>::E
(
    const scalar& rho,
    const scalar& e
) const
{
    return 0.0;//amma_*(b_ - a_)/max(rho, 1.0e-10);
}


template<class Specie>
Foam::scalar Foam::Tait<Specie>::Cv
(
    const scalar& rho,
    const scalar& e
) const
{
    // Check
    return 0;
}


template<class Specie>
Foam::scalar Foam::Tait<Specie>::Cp
(
    const scalar& rho,
    const scalar& e
) const
{
    // Check
    return 0;
}


template<class Specie>
Foam::scalar Foam::Tait<Specie>::H
(
    const scalar& rho,
    const scalar& e
) const
{
    return 0;
}


template<class Specie>
Foam::scalar Foam::Tait<Specie>::CpMCv
(
    const scalar& rho,
    const scalar& e
) const
{
    return this->R();
}


template<class Specie>
Foam::scalar Foam::Tait<Specie>::S
(
    const scalar& p,
    const scalar& rho,
    const scalar& T
) const
{
    // Check
    return -this->R()*log(p/Foam::constant::thermodynamic::Pstd);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::Tait<Specie>::operator+=
(
    const Tait<Specie>& pf
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(pf);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = pf.Y()/this->Y();

        a_ = Y1*a_ + Y2*pf.a_;
        b_ = Y1*b_ + Y2*pf.b_;
        gamma_ = Y1*gamma_ + Y2*pf.gamma_;
    }
}


template<class Specie>
inline void Foam::Tait<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::Tait<Specie> Foam::operator+
(
    const Tait<Specie>& pf1,
    const Tait<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
      + static_cast<const Specie&>(pf2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return Tait<Specie>
        (
            sp,
            pf1.a_,
            pf1.b_,
            pf1.gamma_
        );
    }
    else
    {
        const scalar Y1 = pf1.Y()/sp.Y();
        const scalar Y2 = pf2.Y()/sp.Y();

        return Tait<Specie>
        (
            sp,
            Y1*pf1.a_ + Y2*pf2.a_,
            Y1*pf1.b_ + Y2*pf2.b_,
            Y1*pf1.gamma_ + Y2*pf2.gamma_
        );
    }
}


template<class Specie>
inline Foam::Tait<Specie> Foam::operator*
(
    const scalar s,
    const Tait<Specie>& pf
)
{
    return Tait<Specie>
    (
        s*static_cast<const Specie&>(pf),
        pf.a_,
        pf.b_,
        pf.gamma_
    );
}


template<class Specie>
inline Foam::Tait<Specie> Foam::operator==
(
    const Tait<Specie>& pf1,
    const Tait<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
     == static_cast<const Specie&>(pf2)
    );

    const scalar Y1 = pf1.Y()/sp.Y();
    const scalar Y2 = pf2.Y()/sp.Y();

    return Tait<Specie>
    (
        sp,
        Y2*pf2.a_    - Y1*pf1.a_,
        Y2*pf2.b_  - Y1*pf1.b_,
        Y2*pf2.gamma_ - Y1*pf1.gamma_
    );
}


// ************************************************************************* //

