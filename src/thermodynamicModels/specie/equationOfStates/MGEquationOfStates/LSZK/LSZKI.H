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
inline Foam::LSZK<Specie>::LSZK
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
inline Foam::LSZK<Specie>::LSZK
(
    const word& name,
    const LSZK<Specie>& pf
)
:
    Specie(name, pf),
    a_(pf.a_),
    b_(pf.b_),
    gamma_(pf.gamma_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::LSZK<Specie>>
Foam::LSZK<Specie>::clone() const
{
    return autoPtr<LSZK<Specie>>
    (
        new LSZK<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::LSZK<Specie>>
Foam::LSZK<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<LSZK<Specie>>
    (
        new LSZK<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
Foam::scalar Foam::LSZK<Specie>::Pi
(
    const scalar rho,
    const scalar e
) const
{
    return -a_*pow(rho, b_);
}


template<class Specie>
Foam::scalar Foam::LSZK<Specie>::delta
(
    const scalar p,
    const scalar rho,
    const scalar e
) const
{
    return -a_*b_*pow(rho, b_ - 1.0)/(gamma_ - 1.0);
}


template<class Specie>
Foam::scalar Foam::LSZK<Specie>::GammaMG
(
    const scalar rho,
    const scalar e
) const
{
    return gamma_;
}


template<class Specie>
Foam::scalar Foam::LSZK<Specie>::dpdv
(
    const scalar rho,
    const scalar e
) const
{
    return -sqr(rho)*(gamma_ - 1.0)*e - a_*b_*pow(rho, b_ + 1.0);
}


template<class Specie>
Foam::scalar Foam::LSZK<Specie>::dpde
(
    const scalar rho,
    const scalar e
) const
{
    return (gamma_ - 1.0)*rho;
}


template<class Specie>
Foam::scalar Foam::LSZK<Specie>::E
(
    const scalar rho,
    const scalar e
) const
{
    return a_*pow(rho, b_ - 1.0)/(b_ - 1.0);
}


template<class Specie>
Foam::scalar Foam::LSZK<Specie>::Cv
(
    const scalar rho,
    const scalar e
) const
{
    // Check
    return 0;
}


template<class Specie>
Foam::scalar Foam::LSZK<Specie>::Cp
(
    const scalar rho,
    const scalar e
) const
{
    // Check
    return 0;
}


template<class Specie>
Foam::scalar Foam::LSZK<Specie>::H
(
    const scalar rho,
    const scalar e
) const
{
    return a_*log(rho);
}


template<class Specie>
Foam::scalar Foam::LSZK<Specie>::CpMCv
(
    const scalar rho,
    const scalar e
) const
{
    return this->R();
}


template<class Specie>
Foam::scalar Foam::LSZK<Specie>::S
(
    const scalar p,
    const scalar rho,
    const scalar T
) const
{
    // check
    return -this->R()*log(p/Foam::constant::thermodynamic::Pstd);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::LSZK<Specie>::operator+=
(
    const LSZK<Specie>& pf
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
inline void Foam::LSZK<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::LSZK<Specie> Foam::operator+
(
    const LSZK<Specie>& pf1,
    const LSZK<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
      + static_cast<const Specie&>(pf2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return LSZK<Specie>
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

        return LSZK<Specie>
        (
            sp,
            Y1*pf1.a_ + Y2*pf2.a_,
            Y1*pf1.b_ + Y2*pf2.b_,
            Y1*pf1.gamma_ + Y2*pf2.gamma_
        );
    }
}


template<class Specie>
inline Foam::LSZK<Specie> Foam::operator*
(
    const scalar s,
    const LSZK<Specie>& pf
)
{
    return LSZK<Specie>
    (
        s*static_cast<const Specie&>(pf),
        pf.a_,
        pf.b_,
        pf.gamma_
    );
}


template<class Specie>
inline Foam::LSZK<Specie> Foam::operator==
(
    const LSZK<Specie>& pf1,
    const LSZK<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
     == static_cast<const Specie&>(pf2)
    );

    const scalar Y1 = pf1.Y()/sp.Y();
    const scalar Y2 = pf2.Y()/sp.Y();

    return LSZK<Specie>
    (
        sp,
        Y2*pf2.a_    - Y1*pf1.a_,
        Y2*pf2.b_  - Y1*pf1.b_,
        Y2*pf2.gamma_ - Y1*pf1.gamma_
    );
}


// ************************************************************************* //
