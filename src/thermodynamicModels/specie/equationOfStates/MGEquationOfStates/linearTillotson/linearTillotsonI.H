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
inline Foam::linearTillotson<Specie>::linearTillotson
(
    const Specie& sp,
    const scalar p0,
    const scalar rho0,
    const scalar e0,
    const scalar omega,
    const scalar A,
    const scalar B,
    const scalar C,
    const scalar pCav
)
:
    Specie(sp),
    p0_(p0),
    rho0_(rho0),
    e0_(e0),
    omega_(omega),
    A_(A),
    B_(B),
    C_(C),
    pCav_(pCav)
{}


template<class Specie>
Foam::scalar Foam::linearTillotson<Specie>::Pi
(
    const scalar rho,
    const scalar e
) const
{
    scalar mu(rho/rho0_ - 1.0);
    return
        omega_*rho*e0_
      - p0_
      - A_*mu
      - B_*sqr(mu)
      - C_*pow3(mu);
}


template<class Specie>
Foam::scalar Foam::linearTillotson<Specie>::delta
(
    const scalar p,
    const scalar rho,
    const scalar e
) const
{
    scalar mu(rho/rho0_ - 1.0);
    return
        e0_
      - (
            A_
          + 2.0*B_*mu
          + 3.0*C_*sqr(mu)
        )/(rho0_*omega_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::linearTillotson<Specie>::linearTillotson
(
    const word& name,
    const linearTillotson<Specie>& pf
)
:
    Specie(name, pf),
    p0_(pf.p0_),
    rho0_(pf.rho0_),
    e0_(pf.e0_),
    omega_(pf.omega_),
    A_(pf.A_),
    B_(pf.B_),
    C_(pf.C_),
    pCav_(pf.pCav_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::linearTillotson<Specie>>
Foam::linearTillotson<Specie>::clone() const
{
    return autoPtr<linearTillotson<Specie>>
    (
        new linearTillotson<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::linearTillotson<Specie>>
Foam::linearTillotson<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<linearTillotson<Specie>>
    (
        new linearTillotson<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
Foam::scalar Foam::linearTillotson<Specie>::p
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const bool limit
) const
{
    scalar mu(rho/rho0_ - 1.0);
    scalar pressure =
        p0_
      + omega_*rho*(e - e0_)
      + A_*mu
      + B_*sqr(mu)
      + C_*pow3(mu);

    return limit ? max(pressure, pCav_) : pressure;
}


template<class Specie>
Foam::scalar Foam::linearTillotson<Specie>::Gamma
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
Foam::scalar Foam::linearTillotson<Specie>::cSqr
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
Foam::scalar Foam::linearTillotson<Specie>::dpdv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    scalar mu(rho/rho0_ - 1.0);
    // dRho/dv*dp/dRho
    return
      - sqr(rho)*omega_*(e - e0_)
      - sqr(rho)/rho0_*(A_ + 2.0*B_*mu + 3.0*C_*sqr(mu));
}


template<class Specie>
Foam::scalar Foam::linearTillotson<Specie>::dpde
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return omega_*rho;
}


template<class Specie>
Foam::scalar Foam::linearTillotson<Specie>::dpdT
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
Foam::scalar Foam::linearTillotson<Specie>::E
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return EcTable_.lookup(rho) - e0_;
}


template<class Specie>
Foam::scalar Foam::linearTillotson<Specie>::Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::linearTillotson<Specie>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::linearTillotson<Specie>::CpMCv
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
Foam::scalar Foam::linearTillotson<Specie>::H
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::linearTillotson<Specie>::S
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
inline void Foam::linearTillotson<Specie>::operator+=
(
    const linearTillotson<Specie>& pf
)
{
    NotImplemented;
    scalar Y1 = this->Y();
    Specie::operator+=(pf);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = pf.Y()/this->Y();

        p0_ = Y1*p0_ + Y2*pf.p0_;
        rho0_ = Y1*rho0_ + Y2*pf.rho0_;
        e0_ = Y1*e0_ + Y2*pf.e0_;
        omega_ = Y1*omega_ + Y2*pf.omega_;
        A_ = Y1*A_ + Y2*pf.A_;
        B_ = Y1*B_ + Y2*pf.B_;
        C_ = Y1*C_ + Y2*pf.C_;
        pCav_ = Y1*pCav_ + Y2*pf.pCav_;
    }
}


template<class Specie>
inline void Foam::linearTillotson<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::linearTillotson<Specie> Foam::operator+
(
    const linearTillotson<Specie>& pf1,
    const linearTillotson<Specie>& pf2
)
{
    NotImplemented;
    Specie sp
    (
        static_cast<const Specie&>(pf1)
      + static_cast<const Specie&>(pf2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return linearTillotson<Specie>
        (
            sp,
            pf1.p0_,
            pf1.rho0_,
            pf1.e0_,
            pf1.omega_,
            pf1.A_,
            pf1.B_,
            pf1.C_,
            pf1.pCav_
        );
    }
    else
    {
        const scalar Y1 = pf1.Y()/sp.Y();
        const scalar Y2 = pf2.Y()/sp.Y();

        return linearTillotson<Specie>
        (
            sp,
            Y1*pf1.p0_ + Y2*pf2.p0_,
            Y1*pf1.rho0_ + Y2*pf2.rho0_,
            Y1*pf1.e0_ + Y2*pf2.e0_,
            Y1*pf1.omega_ + Y2*pf2.omega_,
            Y1*pf1.A_ + Y2*pf2.A_,
            Y1*pf1.B_ + Y2*pf2.B_,
            Y1*pf1.C_ + Y2*pf2.C_,
            Y1*pf1.pCav_ + Y2*pf2.pCav_
        );
    }
}


template<class Specie>
inline Foam::linearTillotson<Specie> Foam::operator*
(
    const scalar s,
    const linearTillotson<Specie>& pf
)
{
    return linearTillotson<Specie>
    (
        s*static_cast<const Specie&>(pf),
        pf.p0_,
        pf.rho0_,
        pf.e0_,
        pf.omega_,
        pf.A_,
        pf.B_,
        pf.C_,
        pf.pCav_
    );
}


template<class Specie>
inline Foam::linearTillotson<Specie> Foam::operator==
(
    const linearTillotson<Specie>& pf1,
    const linearTillotson<Specie>& pf2
)
{
    NotImplemented;
    Specie sp
    (
        static_cast<const Specie&>(pf1)
     == static_cast<const Specie&>(pf2)
    );

    const scalar Y1 = pf1.Y()/sp.Y();
    const scalar Y2 = pf2.Y()/sp.Y();

    return linearTillotson<Specie>
    (
        sp,
        Y2*pf2.p0_    - Y1*pf1.p0_,
        Y2*pf2.rho0_    - Y1*pf1.rho0_,
        Y2*pf2.e0_    - Y1*pf1.e0_,
        Y2*pf2.omega_  - Y1*pf1.omega_,
        Y2*pf2.A_  - Y1*pf1.A_,
        Y2*pf2.B_ - Y1*pf1.B_,
        Y2*pf2.C_     - Y1*pf1.C_,
        Y2*pf2.pCav_  - Y1*pf1.pCav_
    );
}


// ************************************************************************* //
