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
inline Foam::AbelNobel<Specie>::AbelNobel
(
    const Specie& sp,
    const scalar b
)
:
    Specie(sp),
    b_(b)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::AbelNobel<Specie>::AbelNobel
(
    const word& name,
    const AbelNobel<Specie>& pf
)
:
    Specie(name, pf),
    b_(pf.b_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::AbelNobel<Specie>>
Foam::AbelNobel<Specie>::clone() const
{
    return autoPtr<AbelNobel<Specie>>
    (
        new AbelNobel<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::AbelNobel<Specie>>
Foam::AbelNobel<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<AbelNobel<Specie>>
    (
        new AbelNobel<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
Foam::scalar Foam::AbelNobel<Specie>::p
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
Foam::scalar Foam::AbelNobel<Specie>::pRhoT
(
    const scalar rho,
    const scalar T
) const
{
    return rho*this->R()*T/(1.0 - b_*rho);
}


template<class Specie>
Foam::scalar Foam::AbelNobel<Specie>::rho
(
    const scalar p,
    const scalar T
) const
{
    return p/(this->R()*T + p*b_);
}


template<class Specie>
Foam::scalar Foam::AbelNobel<Specie>::dpdv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return -this->R()*T*sqr(rho/(b_*rho - 1.0));
}


template<class Specie>
Foam::scalar Foam::AbelNobel<Specie>::dpde
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
Foam::scalar Foam::AbelNobel<Specie>::dpdT
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return rho*this->R()/(1.0 - b_*rho);
}


template<class Specie>
Foam::scalar Foam::AbelNobel<Specie>::Gamma
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar cv
) const
{
    return max(dpdT(rho, e, T)/(cv*max(rho, 1e-10)), small) + 1.0;
}


template<class Specie>
Foam::scalar Foam::AbelNobel<Specie>::cSqr
(
    const scalar p,
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar cv
) const
{
    return
        (sqr(dpdT(rho, e, T))*T/cv - dpdv(rho, e, T))
       /sqr(max(rho, 1e-10));
}


template<class Specie>
Foam::scalar Foam::AbelNobel<Specie>::E
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::AbelNobel<Specie>::Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::AbelNobel<Specie>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::AbelNobel<Specie>::H
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::AbelNobel<Specie>::CpMCv
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
Foam::scalar Foam::AbelNobel<Specie>::S
(
    const scalar p,
    const scalar T
) const
{
    return 0.0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::AbelNobel<Specie>::operator+=
(
    const AbelNobel<Specie>& pf
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(pf);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = pf.Y()/this->Y();

        b_ = Y1*b_ + Y2*pf.b_;
    }
}


template<class Specie>
inline void Foam::AbelNobel<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::AbelNobel<Specie> Foam::operator+
(
    const AbelNobel<Specie>& pf1,
    const AbelNobel<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
      + static_cast<const Specie&>(pf2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return AbelNobel<Specie>
        (
            sp,
            pf1.b_
        );
    }
    else
    {
        const scalar Y1 = pf1.Y()/sp.Y();
        const scalar Y2 = pf2.Y()/sp.Y();

        return AbelNobel<Specie>
        (
            sp,
            Y1*pf1.b_ + Y2*pf2.b_
        );
    }
}


template<class Specie>
inline Foam::AbelNobel<Specie> Foam::operator*
(
    const scalar s,
    const AbelNobel<Specie>& pf
)
{
    return AbelNobel<Specie>
    (
        s*static_cast<const Specie&>(pf),
        pf.b_
    );
}


template<class Specie>
inline Foam::AbelNobel<Specie> Foam::operator==
(
    const AbelNobel<Specie>& pf1,
    const AbelNobel<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
     == static_cast<const Specie&>(pf2)
    );

    const scalar Y1 = pf1.Y()/sp.Y();
    const scalar Y2 = pf2.Y()/sp.Y();

    return AbelNobel<Specie>
    (
        sp,
        Y2*pf2.b_ - Y1*pf1.b_
    );
}


// ************************************************************************* //
