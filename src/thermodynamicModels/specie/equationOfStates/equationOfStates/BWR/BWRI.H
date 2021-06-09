/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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
inline Foam::BWR<Specie>::BWR
(
    const Specie& sp,
    const scalar A0,
    const scalar B0,
    const scalar C0,
    const scalar a,
    const scalar b,
    const scalar c,
    const scalar alpha,
    const scalar rhoc,
    const scalar gamma
)
:
    Specie(sp),
    A0_(A0),
    B0_(B0),
    C0_(C0),
    a_(a),
    b_(b),
    c_(c),
    alpha_(alpha),
    rhoc_(rhoc),
    gamma_(gamma)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::BWR<Specie>::BWR
(
    const word& name,
    const BWR<Specie>& pf
)
:
    Specie(name, pf),
    A0_(pf.A0_),
    B0_(pf.B0_),
    C0_(pf.C0_),
    a_(pf.a_),
    b_(pf.b_),
    c_(pf.c_),
    alpha_(pf.alpha_),
    rhoc_(pf.rhoc_),
    gamma_(pf.gamma)
{}


template<class Specie>
inline Foam::autoPtr<Foam::BWR<Specie>>
Foam::BWR<Specie>::clone() const
{
    return autoPtr<BWR<Specie>>
    (
        new BWR<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::BWR<Specie>>
Foam::BWR<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<BWR<Specie>>
    (
        new BWR<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
Foam::scalar Foam::BWR<Specie>::pRhoT
(
    const scalar rho,
    const scalar T
) const
{
    return
      - this->R()*T*(sqr(rho) + 2.0*pow3(rho)*B0_ + 3.0*pow4(rho)*b_)
      + 2.0*A0_*pow3(rho)
      + 3.0*a_*pow4(rho)
      - 6.0*a_*alpha_*pow(rho, 7)
      + 2.0*C0_*pow3(rho)/sqr(T)
      + c_*pow4(rho)/sqr(T)*exp(-gamma_*sqr(rho))
       *(2.0*sqr(gamma_*sqr(rho)) - 3.0*gamma_*sqr(rho) - 3.0);
}


template<class Specie>
Foam::scalar Foam::BWR<Specie>::dpdv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return
      - this->R()*T
       *(sqr(rho) + 2.0*pow3(rho)*B0_ + 3.0*pow4(rho)*b_)
      + 2.0*A0_*pow3(rho)
      + 3.0*a_*pow4(rho)
      - 6.0*a_*alpha_*pow(rho, 7)
      + 2.0*C0_*pow3(rho)/sqr(T)
      + c_*pow4(rho)/sqr(T)*exp(-gamma_*sqr(rho))
       *(2.0*sqr(gamma_*sqr(rho)) - 3.0*gamma_*sqr(rho) - 3.0);
}


template<class Specie>
Foam::scalar Foam::BWR<Specie>::dpdT
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return
        this->R()
       *(
            rho
          + sqr(rho)*B0_
          + pow3(rho)*b_
        )
      + 2.0*C0_*sqr(rho)/pow3(T)
      - 2.0*pow3(rho/T)*c_*(1.0 + gamma_*sqr(rho))*exp(-gamma_*sqr(rho));
}


template<class Specie>
Foam::scalar Foam::BWR<Specie>::E
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return
      - A0_*rho
      - 3.0*C0_*rho/sqr(T)
      + a_*alpha_*pow5(rho)/5.0
      - a_*sqr(rho)/2.0
      - 3.0*c_/sqr(T)*exp(-gamma_*sqr(rho))*(sqr(rho)/2.0 + 1.0/gamma_);
}


template<class Specie>
Foam::scalar Foam::BWR<Specie>::Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return
        (6.0*C0_*gamma_*rho + c_*(gamma_*sqr(rho) + 2.0)*exp(-gamma_*sqr(rho)))
       /(pow4(T)*gamma_);
}


template<class Specie>
Foam::scalar Foam::BWR<Specie>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Cv(rho, T) - T*sqr(dpdT(rho, T))/dpdv(rho, T) - this->R();
}


template<class Specie>
Foam::scalar Foam::BWR<Specie>::H
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return
        A0_*rho
      + 3.0*C0_*rho/sqr(T)
      - this->R()*T/2.0*(b_*sqr(rho) - 1)
      - 3.0*a_*alpha_*pow5(rho)
      + 1.5*a_*sqr(rho)
      + c_*sqr(rho)*exp(-gamma_*sqr(rho))/sqr(T)
       *(sqr(gamma_)*pow4(rho) - 3.5*(gamma_*sqr(rho) + 1.0));
}


template<class Specie>
Foam::scalar Foam::BWR<Specie>::CpMCv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return
        2.0*C0_*rho/pow3(T)
      + this->R()*(B0_*rho + b_*sqr(rho) + 1.0)
      - 2.0*c_*sqr(rho)*(gamma_*sqr(rho) + 1.0)*exp(-gamma_*sqr(rho))/pow3(T);
}


template<class Specie>
Foam::scalar Foam::BWR<Specie>::S
(
    const scalar p,
    const scalar rho,
    const scalar T
) const
{
    return 0.0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::BWR<Specie>::operator+=
(
    const BWR<Specie>& pf
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(pf);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = pf.Y()/this->Y();

        A0_ = Y1*A0_ + Y2*pf.A0_;
        B0_ = Y1*B0_ + Y2*pf.B0_;
        C0_ = Y1*C0_ + Y2*pf.C0_;
        a_ = Y1*a_ + Y2*pf.a_;
        b_ = Y1*b_ + Y2*pf.b_;
        c_ = Y1*c_ + Y2*pf.c_;
        alpha_ = Y1*alpha_ + Y2*pf.alpha_;
        rhoc_ = Y1*rhoc_ + Y2*pf.rhoc_;
        gamma_ = Y1*gamma_ + Y2*pf.gamma_;
    }
}


template<class Specie>
inline void Foam::BWR<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::BWR<Specie> Foam::operator+
(
    const BWR<Specie>& pf1,
    const BWR<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
      + static_cast<const Specie&>(pf2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return BWR<Specie>
        (
            sp,
            pf1.A0_,
            pf1.B0_,
            pf1.C0_,
            pf1.a_,
            pf1.b_,
            pf1.c_,
            pf1.alpha_,
            pf1.rhoc_,
            pf1.gamma_
        );
    }
    else
    {
        const scalar Y1 = pf1.Y()/sp.Y();
        const scalar Y2 = pf2.Y()/sp.Y();

        return BWR<Specie>
        (
            sp,
            Y1*pf1.A0_ + Y2*pf2.A0_,
            Y1*pf1.B0_ + Y2*pf2.B0_,
            Y1*pf1.C0_ + Y2*pf2.C0_,
            Y1*pf1.a_ + Y2*pf2.a_,
            Y1*pf1.b_ + Y2*pf2.b_,
            Y1*pf1.c_ + Y2*pf2.c_,
            Y1*pf1.alpha_ + Y2*pf2.alpha_,
            Y1*pf1.rhoc_ + Y2*pf2.rhoc_,
            Y1*pf1.gamma_ + Y2*pf2.gamma_
        );
    }
}


template<class Specie>
inline Foam::BWR<Specie> Foam::operator*
(
    const scalar s,
    const BWR<Specie>& pf
)
{
    return BWR<Specie>
    (
        s*static_cast<const Specie&>(pf),
        pf.A0_,
        pf.B0_,
        pf.C0_,
        pf.a_,
        pf.b_,
        pf.c_,
        pf.alpha_,
        pf.rhoc_,
        pf.gamma_
    );
}


template<class Specie>
inline Foam::BWR<Specie> Foam::operator==
(
    const BWR<Specie>& pf1,
    const BWR<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
     == static_cast<const Specie&>(pf2)
    );

    const scalar Y1 = pf1.Y()/sp.Y();
    const scalar Y2 = pf2.Y()/sp.Y();

    return BWR<Specie>
    (
        sp,
        Y2*pf2.A0_ - Y1*pf1.A0_,
        Y2*pf2.B0_ - Y1*pf1.B0_,
        Y2*pf2.C0_ - Y1*pf1.C0_,
        Y2*pf2.a_ - Y1*pf1.a_,
        Y2*pf2.b_ - Y1*pf1.b_,
        Y2*pf2.c_ - Y1*pf1.c_,
        Y2*pf2.alpha_ - Y1*pf1.alpha_,
        Y2*pf2.rhoc_ - Y1*pf1.rhoc_,
        Y2*pf2.gamma_ - Y1*pf1.gamma_
    );
}


// ************************************************************************* //
