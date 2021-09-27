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
inline Foam::CochranChan<Specie>::CochranChan
(
    const Specie& sp,
    const scalar rho0,
    const scalar Gamma0,
    const scalar A,
    const scalar Epsilon1,
    const scalar B,
    const scalar Epsilon2,
    const scalar e0
)
:
    Specie(sp),
    rho0_(rho0),
    Gamma0_(Gamma0),
    A_(A),
    Epsilon1_(Epsilon1),
    B_(B),
    Epsilon2_(Epsilon2),
    e0_(e0)
{}


template<class Specie>
Foam::scalar Foam::CochranChan<Specie>::Pi
(
    const scalar rho,
    const scalar e
) const
{
    scalar rhos(max(rho, 1e-10));
    return
        Gamma0_*rho
       *(
          - A_/((1.0 - Epsilon1_)*rho0_)*(pow(rho0_/rhos, 1.0 - Epsilon1_) - 1.0)
          + B_/((1.0 - Epsilon2_)*rho0_)*(pow(rho0_/rhos, 1.0 - Epsilon2_) - 1.0)
        )
      - A_*pow(rho0_/rhos, -Epsilon1_) + B_*pow(rho0_/rhos, -Epsilon2_)
      - e0_*rho*Gamma0_;
}


template<class Specie>
Foam::scalar Foam::CochranChan<Specie>::delta
(
    const scalar p,
    const scalar rho,
    const scalar e
) const
{
    scalar rhos(max(rho, 1e-10));
    return
      (
           -A_
           *(
               Epsilon1_*pow(rho0_/rhos, -Epsilon1_)
              *(Epsilon1_ - Gamma0_ - 1.0)/rhos
             + Gamma0_/rho0_
            )/(Epsilon1_ - 1.0)
          + B_
           *(
               Epsilon2_*pow(rho0_/rhos, -Epsilon2_)
              *(Epsilon2_ - Gamma0_ - 1.0)/rhos
             + Gamma0_/rho0_
            )/(Epsilon2_ - 1.0)
          - Gamma0_*e0_
        )/Gamma0_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::CochranChan<Specie>::CochranChan
(
    const word& name,
    const CochranChan<Specie>& pf
)
:
    Specie(name, pf),
    rho0_(pf.rho0_),
    Gamma0_(pf.Gamma0_),
    A_(pf.A_),
    Epsilon1_(pf.Epsilon1_),
    B_(pf.B_),
    Epsilon2_(pf.Epsilon2_),
    e0_(pf.e0_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::CochranChan<Specie>>
Foam::CochranChan<Specie>::clone() const
{
    return autoPtr<CochranChan<Specie>>
    (
        new CochranChan<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::CochranChan<Specie>>
Foam::CochranChan<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<CochranChan<Specie>>
    (
        new CochranChan<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
Foam::scalar Foam::CochranChan<Specie>::p
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const bool limit
) const
{
    return
        limit
      ? max(Gamma0_*rho*e - Pi(rho, e), 0.0)
      : Gamma0_*rho*e - Pi(rho, e);
}


template<class Specie>
Foam::scalar Foam::CochranChan<Specie>::Gamma
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar Cv
) const
{
    return Gamma0_ + 1.0;
}


template<class Specie>
Foam::scalar Foam::CochranChan<Specie>::cSqr
(
    const scalar p,
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar Cv
) const
{
    scalar h(((Gamma0_ + 1.0)*p + Pi(rho, e))/(Gamma0_*max(rho, 1e-10)));
    return (h - delta(p, rho, e))*Gamma0_;
}


template<class Specie>
Foam::scalar Foam::CochranChan<Specie>::dpdv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    scalar V(rho0_/max(rho, 1e-10));
    return
      - Gamma0_*(e + e0_)*sqr(rho)
      - rho*Epsilon1_*A_*pow(V, -Epsilon1_)
      - rho*Epsilon2_*B_*pow(V, -Epsilon2_);
}


template<class Specie>
Foam::scalar Foam::CochranChan<Specie>::dpde
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Gamma0_*rho;
}


template<class Specie>
Foam::scalar Foam::CochranChan<Specie>::dpdT
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
Foam::scalar Foam::CochranChan<Specie>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::CochranChan<Specie>::Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::CochranChan<Specie>::CpMCv
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar CpCv
) const
{
    return this->eBased_ ? Gamma0_*CpCv : CpCv*Gamma0_/(Gamma0_ + 1.0);
}


template<class Specie>
Foam::scalar Foam::CochranChan<Specie>::E
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    scalar V(rho0_/max(rho, 1e-10));
    return
        A_/((Epsilon1_ - 1.0)*rho0_)*pow(V, 1.0 - Epsilon1_)
      - B_/((Epsilon2_ - 1.0)*rho0_)*pow(V, 1.0 - Epsilon2_)
      - e0_;
}


template<class Specie>
Foam::scalar Foam::CochranChan<Specie>::H
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    scalar V(rho0_/max(rho, 1e-10));
    return
        (
            Epsilon1_*(Epsilon1_ - Gamma0_ - 1.0)*V*pow(V, -Epsilon1_)
          + (Epsilon1_ - 1.0)*Gamma0_*log(rho)
        )*A_/(rho0_*sqr(Epsilon1_ - 1.0))
      + (
            Epsilon2_*(Epsilon2_ - Gamma0_ - 1.0)*V*pow(V, -Epsilon2_)
          + (Epsilon2_ - 1.0)*Gamma0_*log(rho)
        )*B_/(rho0_*sqr(Epsilon2_ - 1.0));
}


template<class Specie>
Foam::scalar Foam::CochranChan<Specie>::S
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
inline void Foam::CochranChan<Specie>::operator+=
(
    const CochranChan<Specie>& pf
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(pf);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = pf.Y()/this->Y();

        rho0_ = Y1*rho0_ + Y2*pf.rho0_;
        Gamma0_ = Y1*Gamma0_ + Y2*pf.Gamma0_;
        A_ = Y1*A_ + Y2*pf.A_;
        Epsilon1_ = Y1*Epsilon1_ + Y2*pf.Epsilon1_;
        B_ = Y1*B_ + Y2*pf.B_;
        Epsilon2_ = Y1*Epsilon2_ + Y2*pf.Epsilon2_;
        e0_= Y1*e0_ + Y2*pf.e0_;
    }
}


template<class Specie>
inline void Foam::CochranChan<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::CochranChan<Specie> Foam::operator+
(
    const CochranChan<Specie>& pf1,
    const CochranChan<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
      + static_cast<const Specie&>(pf2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return CochranChan<Specie>
        (
            sp,
            pf1.rho0_,
            pf1.Gamma0_,
            pf1.A_,
            pf1.Epsilon1_,
            pf1.B_,
            pf1.Epsilon2_,
            pf1.e0_
        );
    }
    else
    {
        const scalar Y1 = pf1.Y()/sp.Y();
        const scalar Y2 = pf2.Y()/sp.Y();

        return CochranChan<Specie>
        (
            sp,
            Y1*pf1.rho0_ + Y2*pf2.rho0_,
            Y1*pf1.Gamma0_ + Y2*pf2.Gamma0_,
            Y1*pf1.A_ + Y2*pf2.A_,
            Y1*pf1.Epsilon1_ + Y2*pf2.Epsilon1_,
            Y1*pf1.B_ + Y2*pf2.B_,
            Y1*pf1.Epsilon2_ + Y2*pf2.Epsilon2_,
            Y1*pf1.e0_ + Y2*pf2.e0_
        );
    }
}


template<class Specie>
inline Foam::CochranChan<Specie> Foam::operator*
(
    const scalar s,
    const CochranChan<Specie>& pf
)
{
    return CochranChan<Specie>
    (
        s*static_cast<const Specie&>(pf),
        pf.rho0_,
        pf.Gamma0_,
        pf.A_,
        pf.Epsilon1_,
        pf.B_,
        pf.Epsilon2_,
        pf.e0_
    );
}


template<class Specie>
inline Foam::CochranChan<Specie> Foam::operator==
(
    const CochranChan<Specie>& pf1,
    const CochranChan<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
     == static_cast<const Specie&>(pf2)
    );

    const scalar Y1 = pf1.Y()/sp.Y();
    const scalar Y2 = pf2.Y()/sp.Y();

    return CochranChan<Specie>
    (
        sp,
        Y2*pf2.rho0_ - Y1*pf1.rho0_,
        Y2*pf2.Gamma0_ - Y1*pf1.Gamma0_,
        Y2*pf2.A_ - Y1*pf1.A_,
        Y2*pf2.Epsilon1_ - Y1*pf1.Epsilon1_,
        Y2*pf2.B_ - Y1*pf1.B_,
        Y2*pf2.Epsilon2_ - Y1*pf1.Epsilon2_,
        Y2*pf2.e0_ - Y1*pf1.e0_
    );
}


// ************************************************************************* //
