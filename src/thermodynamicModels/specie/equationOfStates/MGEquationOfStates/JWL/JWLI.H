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
inline Foam::JWL<Specie>::JWL
(
    const Specie& sp,
    const scalar rho0,
    const scalar rhoCutOff,
    const scalar omega,
    const scalar gammaIdeal,
    const scalar A,
    const scalar B,
    const scalar R1,
    const scalar R2,
    const scalar e0
)
:
    Specie(sp),
    rho0_(rho0),
    rhoCutOff_(rhoCutOff),
    omega_(omega),
    gammaIdeal_(gammaIdeal),
    A_(A),
    B_(B),
    R1_(R1),
    R2_(R2),
    e0_(e0)
{}


template<class Specie>
Foam::scalar Foam::JWL<Specie>::Pi
(
    const scalar rho,
    const scalar e
) const
{
    if (rho < rhoCutOff_)
    {
        return 0.0;
    }

    scalar V(rho0_/max(rho, 1e-10));
    return
        A_*(omega_/(R1_*V) - 1.0)*exp(-R1_*V)
      + B_*(omega_/(R2_*V) - 1.0)*exp(-R2_*V)
      - omega_*rho*e0_;

}


template<class Specie>
Foam::scalar Foam::JWL<Specie>::delta
(
    const scalar p,
    const scalar rho,
    const scalar e
) const
{
    if (rho < rhoCutOff_)
    {
        return 0.0;
    }

    scalar rhos(max(rho, 1e-10));
    return
        (
            A_*exp(-R1_*rho0_/rhos)
           *(omega_*(1.0/(R1_*rho0_) + 1.0/rhos) - R1_*rho0_/sqr(rhos))
          + B_*exp(-R2_*rho0_/rhos)
           *(omega_*(1.0/(R2_*rho0_) + 1.0/rhos) - R2_*rho0_/sqr(rhos))
          - omega_*e0_
        )/omega_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::JWL<Specie>::JWL
(
    const word& name,
    const JWL<Specie>& pf
)
:
    Specie(name, pf),
    rho0_(pf.rho0_),
    rhoCutOff_(pf.rhoCutOff_),
    omega_(pf.omega_),
    gammaIdeal_(pf.gammaIdeal_),
    A_(pf.A_),
    B_(pf.B_),
    R1_(pf.R1_),
    R2_(pf.R2_),
    e0_(pf.e0_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::JWL<Specie>>
Foam::JWL<Specie>::clone() const
{
    return autoPtr<JWL<Specie>>
    (
        new JWL<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::JWL<Specie>>
Foam::JWL<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<JWL<Specie>>
    (
        new JWL<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
Foam::scalar Foam::JWL<Specie>::p
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const bool limit
) const
{
    return
        limit
      ? max(omega_*rho*e - Pi(rho, e), 0.0)
      : omega_*rho*e - Pi(rho, e);
}


template<class Specie>
Foam::scalar Foam::JWL<Specie>::Gamma
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar Cv
) const
{
    return omega_ + 1.0;
}


template<class Specie>
Foam::scalar Foam::JWL<Specie>::cSqr
(
    const scalar p,
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar Cv
) const
{
    scalar h(((omega_ + 1.0)*p + Pi(rho, e))/(omega_*max(rho, 1e-10)));
    return (h - delta(p, rho, e))*omega_;
}


template<class Specie>
Foam::scalar Foam::JWL<Specie>::dpdv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (rho < rhoCutOff_)
    {
        return (gammaIdeal_ - 1.0)*e*sqr(rho);
    }

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
Foam::scalar Foam::JWL<Specie>::dpde
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (rho < rhoCutOff_)
    {
        return (gammaIdeal_ - 1.0)*rho;
    }
    return omega_*rho;
}


template<class Specie>
Foam::scalar Foam::JWL<Specie>::dpdT
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
Foam::scalar Foam::JWL<Specie>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::JWL<Specie>::Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::JWL<Specie>::CpMCv
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar CpCv
) const
{
    return this->eBased_ ? omega_*CpCv : CpCv*omega_/(omega_ + 1.0);
}


template<class Specie>
Foam::scalar Foam::JWL<Specie>::E
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (rho < rhoCutOff_)
    {
        return 0.0;
    }

    scalar rhos = max(rho, 1e-10);
    scalar V(rho0_/rhos);
    return (A_/R1_*exp(-R1_*V) + B_/R2_*exp(-R2_*V))/rho0_ - e0_;
}


template<class Specie>
Foam::scalar Foam::JWL<Specie>::H
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (rho < rhoCutOff_)
    {
        return 0.0;
    }

    scalar rhos(max(rho, 1e-10));
    return
        A_*(1.0/rhos + 1.0/(R1_*rho0_))*exp(-R1_*rho0_/rhos)
      + B_*(1.0/rhos + 1.0/(R2_*rho0_))*exp(-R2_*rho0_/rhos)
      + e0_*omega_*log(rhos);
}


template<class Specie>
Foam::scalar Foam::JWL<Specie>::S
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
inline void Foam::JWL<Specie>::operator+=
(
    const JWL<Specie>& pf
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(pf);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = pf.Y()/this->Y();

        rho0_ = Y1*rho0_ + Y2*pf.rho0_;
        rhoCutOff_ = Y1*rhoCutOff_ + Y2*pf.rhoCutOff_;
        omega_ = Y1*omega_ + Y2*pf.omega_;
        gammaIdeal_ = Y1*gammaIdeal_ + Y2*pf.gammaIdeal_;
        A_ = Y1*A_ + Y2*pf.A_;
        B_ = Y1*B_ + Y2*pf.B_;
        R1_ = Y1*R1_ + Y2*pf.R1_;
        R2_ = Y1*R2_ + Y2*pf.R2_;
        e0_ = Y1*e0_ + Y2*pf.e0_;
    }
}


template<class Specie>
inline void Foam::JWL<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::JWL<Specie> Foam::operator+
(
    const JWL<Specie>& pf1,
    const JWL<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
      + static_cast<const Specie&>(pf2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return JWL<Specie>
        (
            sp,
            pf1.rho0_,
            pf1.rhoCutOff_,
            pf1.omega_,
            pf1.gammaIdeal_,
            pf1.A_,
            pf1.B_,
            pf1.R1_,
            pf1.R2_,
            pf1.e0_
        );
    }
    else
    {
        const scalar Y1 = pf1.Y()/sp.Y();
        const scalar Y2 = pf2.Y()/sp.Y();

        return JWL<Specie>
        (
            sp,
            Y1*pf1.rho0_ + Y2*pf2.rho0_,
            Y1*pf1.rhoCutOff_ + Y2*pf2.rhoCutOff_,
            Y1*pf1.omega_ + Y2*pf2.omega_,
            Y1*pf1.gammaIdeal_ + Y2*pf2.gammaIdeal_,
            Y1*pf1.A_ + Y2*pf2.A_,
            Y1*pf1.B_ + Y2*pf2.B_,
            Y1*pf1.R1_ + Y2*pf2.R1_,
            Y1*pf1.R2_ + Y2*pf2.R2_,
            Y1*pf1.e0_ + Y2*pf2.e0_
        );
    }
}


template<class Specie>
inline Foam::JWL<Specie> Foam::operator*
(
    const scalar s,
    const JWL<Specie>& pf
)
{
    return JWL<Specie>
    (
        s*static_cast<const Specie&>(pf),
        pf.rho0_,
        pf.rhoCutOff_,
        pf.omega_,
        pf.gammaIdeal_,
        pf.A_,
        pf.B_,
        pf.R1_,
        pf.R2_,
        pf.e0_
    );
}


template<class Specie>
inline Foam::JWL<Specie> Foam::operator==
(
    const JWL<Specie>& pf1,
    const JWL<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
     == static_cast<const Specie&>(pf2)
    );

    const scalar Y1 = pf1.Y()/sp.Y();
    const scalar Y2 = pf2.Y()/sp.Y();

    return JWL<Specie>
    (
        sp,
        Y2*pf2.rho0_ - Y1*pf1.rho0_,
        Y2*pf2.rhoCutOff_ - Y1*pf1.rhoCutOff_,
        Y2*pf2.omega_ - Y1*pf1.omega_,
        Y2*pf2.gammaIdeal_ - Y1*pf1.gammaIdeal_,
        Y2*pf2.A_ - Y1*pf1.A_,
        Y2*pf2.B_ - Y1*pf1.B_,
        Y2*pf2.R1_ - Y1*pf1.R1_,
        Y2*pf2.R2_ - Y1*pf1.R2_,
        Y2*pf2.e0_ - Y1*pf1.e0_
    );
}



// ************************************************************************* //
