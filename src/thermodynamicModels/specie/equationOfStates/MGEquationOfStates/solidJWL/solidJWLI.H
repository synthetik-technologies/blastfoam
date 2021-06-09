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
inline Foam::solidJWL<Specie>::solidJWL
(
    const Specie& sp,
    const scalar rho0,
    const scalar omega,
    const scalar A,
    const scalar B,
    const scalar R1,
    const scalar R2,
    const scalar pRef,
    const scalar e0
)
:
    Specie(sp),
    rho0_(rho0),
    omega_(omega),
    A_(A),
    B_(B),
    R1_(R1),
    R2_(R2),
    pRef_(pRef),
    e0_(e0)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::solidJWL<Specie>::solidJWL
(
    const word& name,
    const solidJWL<Specie>& pf
)
:
    Specie(name, pf),
    rho0_(pf.rho0_),
    omega_(pf.omega_),
    A_(pf.A_),
    B_(pf.B_),
    R1_(pf.R1_),
    R2_(pf.R2_),
    pRef_(pf.pRef_),
    e0_(pf.e0_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::solidJWL<Specie>>
Foam::solidJWL<Specie>::clone() const
{
    return autoPtr<solidJWL<Specie>>
    (
        new solidJWL<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::solidJWL<Specie>>
Foam::solidJWL<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<solidJWL<Specie>>
    (
        new solidJWL<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Specie>
Foam::scalar Foam::solidJWL<Specie>::GammaMG
(
    const scalar rho,
    const scalar e
) const
{
    return omega_ + 1.0;
}


template<class Specie>
Foam::scalar Foam::solidJWL<Specie>::Pi
(
    const scalar rho,
    const scalar e
) const
{
    scalar V(rho0_/max(rho, 1e-10));
    return
        A_*(omega_/(R1_*V) - 1.0)*exp(-R1_*V)
      + B_*(omega_/(R2_*V) - 1.0)*exp(-R2_*V)
      - omega_*rho*e0_
      - pRef_;

}


template<class Specie>
Foam::scalar Foam::solidJWL<Specie>::delta
(
    const scalar p,
    const scalar rho,
    const scalar e
) const
{
    scalar rhos(max(rho, 1e-10));
    return
        (
            A_*exp(-R1_*rho0_/rhos)
           *(omega_*(1.0/(R1_*rho0_) + 1.0/rhos) - R1_*rho0_/sqr(rhos))
          + B_*exp(-R2_*rho0_/rhos)
           *(omega_*(1.0/(R2_*rho0_) + 1.0/rhos) - R2_*rho0_/sqr(rhos))
          + omega_*e0_
        )/omega_;
}


template<class Specie>
Foam::scalar Foam::solidJWL<Specie>::dpdv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    scalar rhos(max(rho, 1e-10));
    return
      - (
            A_*exp(-R1_*rho0_/rhos)
           *(R1_*rho0_ - omega_*rho*(1.0 + rho/R1_/rho0_))
          + B_*exp(-R2_*rho0_/rhos)
           *(R2_*rho0_ - omega_*rho*(1.0 + rho/R2_/rho0_))
          + omega_*(e + e0_)*sqr(rhos)
        );
}


template<class Specie>
Foam::scalar Foam::solidJWL<Specie>::dpde
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return omega_*rho;
}


template<class Specie>
Foam::scalar Foam::solidJWL<Specie>::E
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    scalar rhos = max(rho, 1e-10);
    scalar V(rho0_/rhos);
    return (A_/R1_*exp(-R1_*V) + B_/R2_*exp(-R2_*V))/rho0_ - e0_;
}


template<class Specie>
Foam::scalar Foam::solidJWL<Specie>::H
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    scalar rhos(max(rho, 1e-10));
    return
        A_*(1.0/rhos + 1.0/(R1_*rho0_))*exp(-R1_*rho0_/rhos)
      + B_*(1.0/rhos + 1.0/(R2_*rho0_))*exp(-R2_*rho0_/rhos)
      + e0_*omega_*log(rhos);
}


template<class Specie>
Foam::scalar Foam::solidJWL<Specie>::S
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
inline void Foam::solidJWL<Specie>::operator+=
(
    const solidJWL<Specie>& pf
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
        R1_ = Y1*R1_ + Y2*pf.R1_;
        R2_ = Y1*R2_ + Y2*pf.R2_;
        pRef_ = Y1*pRef_ + Y2*pf.pRef_;
        e0_ = Y1*e0_ + Y2*pf.e0_;
    }
}


template<class Specie>
inline void Foam::solidJWL<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::solidJWL<Specie> Foam::operator+
(
    const solidJWL<Specie>& pf1,
    const solidJWL<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
      + static_cast<const Specie&>(pf2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return solidJWL<Specie>
        (
            sp,
            pf1.rho0_,
            pf1.omega_,
            pf1.A_,
            pf1.B_,
            pf1.R1_,
            pf1.R2_,
            pf1.pRef_,
            pf1.e0_
        );
    }
    else
    {
        const scalar Y1 = pf1.Y()/sp.Y();
        const scalar Y2 = pf2.Y()/sp.Y();

        return solidJWL<Specie>
        (
            Y1*pf1.rho0_ + Y2*pf2.rho0_,
            Y1*pf1.omega_ + Y2*pf2.omega_,
            Y1*pf1.A_ + Y2*pf2.A_,
            Y1*pf1.B_ + Y2*pf2.B_,
            Y1*pf1.R1_ + Y2*pf2.R1_,
            Y1*pf1.R2_ + Y2*pf2.R2_,
            Y1*pf1.pref_ + Y2*pf2.pRef_,
            Y1*pf1.e0_ + Y2*pf2.e0_
        );
    }
}


template<class Specie>
inline Foam::solidJWL<Specie> Foam::operator*
(
    const scalar s,
    const solidJWL<Specie>& pf
)
{
    return solidJWL<Specie>
    (
        s*static_cast<const Specie&>(pf),
        pf.rho0_,
        pf.omega_,
        pf.A_,
        pf.B_,
        pf.R1_,
        pf.R2_,
        pf.pRef_,
        pf.e0_
    );
}


template<class Specie>
inline Foam::solidJWL<Specie> Foam::operator==
(
    const solidJWL<Specie>& pf1,
    const solidJWL<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
     == static_cast<const Specie&>(pf2)
    );

    const scalar Y1 = pf1.Y()/sp.Y();
    const scalar Y2 = pf2.Y()/sp.Y();

    return solidJWL<Specie>
    (
        sp,
        Y2*pf2.rho0_ - Y1*pf1.rho0_,
        Y2*pf2.omega_ - Y1*pf1.omega_,
        Y2*pf2.A_ - Y1*pf1.A_,
        Y2*pf2.B_ - Y1*pf1.B_,
        Y2*pf2.R1_ - Y1*pf1.R1_,
        Y2*pf2.R2_ - Y1*pf1.R2_,
        Y2*pf2.pRef_ - Y1*pf1.pRef_,
        Y2*pf2.e0_ - Y1*pf1.e0_
    );
}


// ************************************************************************* //
