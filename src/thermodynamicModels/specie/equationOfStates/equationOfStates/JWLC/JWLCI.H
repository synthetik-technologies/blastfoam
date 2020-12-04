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
inline Foam::JWLC<Specie>::JWLC
(
    const Specie& sp,
    const scalar rho0,
    const scalar omega,
    const scalar A,
    const scalar B,
    const scalar C,
    const scalar R1,
    const scalar R2
)
:
    Specie(sp),
    rho0_(rho0),
    omega_(omega),
    A_(A),
    B_(B),
    C_(C),
    R1_(R1),
    R2_(R2)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::JWLC<Specie>::JWLC
(
    const word& name,
    const JWLC<Specie>& pf
)
:
    Specie(name, pf),
    rho0_(pf.rho0_),
    omega_(pf.omega_),
    A_(pf.A_),
    B_(pf.B_),
    C_(pf.C_),
    R1_(pf.R1_),
    R2_(pf.R2_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::JWLC<Specie>>
Foam::JWLC<Specie>::clone() const
{
    return autoPtr<JWLC<Specie>>
    (
        new JWLC<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::JWLC<Specie>>
Foam::JWLC<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<JWLC<Specie>>
    (
        new JWLC<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
Foam::scalar Foam::JWLC<Specie>::Gamma
(
    const scalar& rho,
    const scalar& T,
    const scalar& cv
) const
{
    return omega_ + 1.0;
}


template<class Specie>
Foam::scalar Foam::JWLC<Specie>::p
(
    const scalar& rho,
    const scalar& T
) const
{
    scalar V(rho0_/max(rho, 1e-10));
    return A_*exp(-R1_*V) + B_*exp(-R2_*V) + C_/pow(V, omega_ + 1.0);
}


template<class Specie>
Foam::scalar Foam::JWLC<Specie>::cSqr
(
    const scalar& p,
    const scalar& rho,
    const scalar& T,
    const scalar& cv
) const
{
    scalar V(rho0_/max(rho, 1e-10));
    return
        sqr(V)/rho0_
       *(
            A_*R1_*exp(-R1_*V)
            + B_*R2_*exp(-R2_*V)
            + C_*(1.0 + omega_)/pow(V, 2.0 + omega_)
        );
}


template<class Specie>
Foam::scalar Foam::JWLC<Specie>::dpdv
(
    const scalar& rho,
    const scalar& T
) const
{
    scalar V(rho0_/max(rho, 1e-10));
    return
        rho0_
       *(
            A_*R1_*exp(-R1_*V)
          + B_*R2_*exp(-R2_*V)
          + C_*(1.0 + omega_)/pow(V, 2.0 + omega_)
        );
}


template<class Specie>
Foam::scalar Foam::JWLC<Specie>::dpdT
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::JWLC<Specie>::E
(
    const scalar& rho,
    const scalar& T
) const
{
    scalar V(rho0_/max(rho, 1e-10));
    return
        (
            A_/R1_*exp(-R1_*V)
          + B_/R2_*exp(-R2_*V)
          + C_/(omega_*pow(V, omega_))
        )/rho0_;
}


template<class Specie>
Foam::scalar Foam::JWLC<Specie>::Cv
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0;
}


template<class Specie>
Foam::scalar Foam::JWLC<Specie>::Cp
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::JWLC<Specie>::H
(
    const scalar& rho,
    const scalar& T
) const
{
    scalar rhos(max(rho, 1e-10));
    scalar V(rho0_/rhos);
    return
      - (
            A_*R1_*rho0_*exp(-R1_*V)
          + B_*R2_*rho0_*exp(-R2_*V)
        )/(2.0*sqr(rhos))
      - C_*(omega_ + 1.0)*pow(V, -omega_)/(2.0*rho0_);
}


template<class Specie>
Foam::scalar Foam::JWLC<Specie>::CpMCv
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::JWLC<Specie>::S
(
    const scalar& p,
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::JWLC<Specie>::operator+=
(
    const JWLC<Specie>& pf
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(pf);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = pf.Y()/this->Y();

        rho0_ = Y1*rho0_ + Y2*pf.rho0_;
        omega_ = Y1*omega_ + Y2*pf.omega_;
        A_ = Y1*A_ + Y2*pf.A_;
        B_ = Y1*B_ + Y2*pf.B_;
        C_ = Y1*C_ + Y2*pf.C_;
        R1_ = Y1*R1_ + Y2*pf.R1_;
        R2_ = Y1*R2_ + Y2*pf.R2_;
    }
}


template<class Specie>
inline void Foam::JWLC<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::JWLC<Specie> Foam::operator+
(
    const JWLC<Specie>& pf1,
    const JWLC<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
      + static_cast<const Specie&>(pf2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return JWLC<Specie>
        (
            sp,
            pf1.rho0_,
            pf1.omega_,
            pf1.A_,
            pf1.B_,
            pf1.C_,
            pf1.R1_,
            pf1.R2_
        );
    }
    else
    {
        const scalar Y1 = pf1.Y()/sp.Y();
        const scalar Y2 = pf2.Y()/sp.Y();

        return JWLC<Specie>
        (
            sp,
            Y1*pf1.rho0_ + Y2*pf2.rho0_,
            Y1*pf1.omega_ + Y2*pf2.omega_,
            Y1*pf1.A_ + Y2*pf2.A_,
            Y1*pf1.B_ + Y2*pf2.B_,
            Y1*pf1.C_ + Y2*pf2.C_,
            Y1*pf1.R1_ + Y2*pf2.R1_,
            Y1*pf1.R2_ + Y2*pf2.R2_
        );
    }
}


template<class Specie>
inline Foam::JWLC<Specie> Foam::operator*
(
    const scalar s,
    const JWLC<Specie>& pf
)
{
    return JWLC<Specie>
    (
        s*static_cast<const Specie&>(pf),
        pf.rho0_,
        pf.omega_,
        pf.A_,
        pf.B_,
        pf.C_,
        pf.R1_,
        pf.R2_
    );
}


template<class Specie>
inline Foam::JWLC<Specie> Foam::operator==
(
    const JWLC<Specie>& pf1,
    const JWLC<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
     == static_cast<const Specie&>(pf2)
    );

    const scalar Y1 = pf1.Y()/sp.Y();
    const scalar Y2 = pf2.Y()/sp.Y();

    return JWLC<Specie>
    (
        sp,
        Y2*pf2.rho0_ - Y1*pf1.rho0_,
        Y2*pf2.omega_ - Y1*pf1.omega_,
        Y2*pf2.A_ - Y1*pf1.A_,
        Y2*pf2.B_ - Y1*pf1.B_,
        Y2*pf2.C_ - Y1*pf1.C_,
        Y2*pf2.R1_ - Y1*pf1.R1_,
        Y2*pf2.R2_ - Y1*pf1.R2_
    );
}


// ************************************************************************* //
