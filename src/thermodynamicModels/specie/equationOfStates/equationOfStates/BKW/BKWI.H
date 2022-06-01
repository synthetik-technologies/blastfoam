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
inline Foam::BKW<Specie>::BKW
(
    const Specie& sp,
    const scalar k,
    const scalar kappa,
    const scalar Theta,
    const scalar a,
    const scalar beta
)
:
    Specie(sp),
    k_(k),
    kappa_(kappa),
    Theta_(Theta),
    a_(a),
    beta_(beta)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::BKW<Specie>::BKW
(
    const word& name,
    const BKW<Specie>& pf
)
:
    Specie(name, pf),
    k_(pf.k_),
    kappa_(pf.kappa_),
    Theta_(pf.Theta_),
    a_(pf.a_),
    beta_(pf.beta_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::BKW<Specie>>
Foam::BKW<Specie>::clone() const
{
    return autoPtr<BKW<Specie>>
    (
        new BKW<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::BKW<Specie>>
Foam::BKW<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<BKW<Specie>>
    (
        new BKW<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
Foam::scalar Foam::BKW<Specie>::p
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
Foam::scalar Foam::BKW<Specie>::pRhoT
(
    const scalar rho,
    const scalar T
) const
{
    scalar x = kappa_*k_*rho/(this->W()*pow(T + Theta_, a_));
    return this->R()*T*rho*(1.0 + x*exp(beta_*x));
}


template<class Specie>
Foam::scalar Foam::BKW<Specie>::rho
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class Specie>
Foam::scalar Foam::BKW<Specie>::Gamma
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
Foam::scalar Foam::BKW<Specie>::cSqr
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
Foam::scalar Foam::BKW<Specie>::dpdv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    scalar x = kappa_*k_*rho/(this->W()*pow(T + Theta_, a_));
    scalar dxdv = -kappa_*k_*sqr(rho)/(this->W()*pow(T + Theta_, a_));
    return
        this->R()*T*rho
       *(
           exp(beta_*x)*dxdv*(1.0 + x*beta_)
         - (1.0 + x*exp(beta_*x))*rho
        );
}


template<class Specie>
Foam::scalar Foam::BKW<Specie>::dpde
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
Foam::scalar Foam::BKW<Specie>::dpdT
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    scalar x = kappa_*k_*rho/(this->W()*pow(T + Theta_, a_));
    scalar dxdT = -a_*kappa_*k_*rho/(this->W()*pow(T + Theta_, a_ + 1.0));
    return
        this->R()
       *(
            1.0 + x*exp(beta_*x)
          + T*exp(beta_*x)*dxdT*(1.0 + x*beta_)
        );
}


template<class Specie>
Foam::scalar Foam::BKW<Specie>::E
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    scalar x = kappa_*k_*rho/(this->W()*pow(T + Theta_, a_));
    return this->R()*sqr(T)*a_*x*exp(beta_*x)/(T + Theta_);
}


template<class Specie>
Foam::scalar Foam::BKW<Specie>::Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    scalar x = kappa_*k_*rho/(this->W()*pow(T + Theta_, a_));
    return
        this->R()*T*x*a_*exp(beta_*x)/sqr(T + Theta_)
       *(1.0 + 2.0*Theta_/T - a_*(1.0 + beta_*x));
}


template<class Specie>
Foam::scalar Foam::BKW<Specie>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Cv(rho, e, T) - T*sqr(dpdT(rho, e, T))/dpdv(rho, e, T) - this->R();
}


template<class Specie>
Foam::scalar Foam::BKW<Specie>::H
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    scalar x = kappa_*k_*rho/(this->W()*pow(T + Theta_, a_));
    return this->R()*T/(T + Theta_)*x*(T*(a_ + 1.0) + Theta_)*exp(beta_*x);
}


template<class Specie>
Foam::scalar Foam::BKW<Specie>::CpMCv
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar CpCv
) const
{
    scalar x = kappa_*k_*rho/(this->W()*pow(T + Theta_, a_));
    return
        this->R()
       *(
          - T*a_*x/(T + Theta_)*(beta_*x + 1.0)*exp(beta_*x)
          + x*exp(beta_*x) + 1
        );
}


template<class Specie>
Foam::scalar Foam::BKW<Specie>::S
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return 0.0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::BKW<Specie>::operator+=
(
    const BKW<Specie>& pf
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(pf);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = pf.Y()/this->Y();

        k_ = Y1*k_ + Y2*pf.k_;
        kappa_ = Y1*kappa_ + Y2*pf.kappa_;
        Theta_ = Y1*Theta_ + Y2*pf.Theta_;
        a_ = Y1*a_ + Y2*pf.a_;
        beta_ = Y1*beta_ + Y2*pf.beta_;
    }
}


template<class Specie>
inline void Foam::BKW<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::BKW<Specie> Foam::operator+
(
    const BKW<Specie>& pf1,
    const BKW<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
      + static_cast<const Specie&>(pf2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return BKW<Specie>
        (
            sp,
            pf1.k_,
            pf1.kappa,
            pf1.Theta_,
            pf1.a_,
            pf1.beta_
        );
    }
    else
    {
        const scalar Y1 = pf1.Y()/sp.Y();
        const scalar Y2 = pf2.Y()/sp.Y();

        return BKW<Specie>
        (
            sp,
            Y1*pf1.k_ + Y2*pf2.k_,
            Y1*pf1.kappa_ + Y2*pf2.kappa,
            Y1*pf1.Theta_ + Y2*pf2.Theta_,
            Y1*pf1.a_ + Y2*pf2.a_,
            Y1*pf1.beta_ + Y2*pf2.beta_
        );
    }
}


template<class Specie>
inline Foam::BKW<Specie> Foam::operator*
(
    const scalar s,
    const BKW<Specie>& pf
)
{
    return BKW<Specie>
    (
        s*static_cast<const Specie&>(pf),
        pf.k_,
        pf.kappa,
        pf.Theta_,
        pf.a_,
        pf.beta_
    );
}


template<class Specie>
inline Foam::BKW<Specie> Foam::operator==
(
    const BKW<Specie>& pf1,
    const BKW<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
     == static_cast<const Specie&>(pf2)
    );

    const scalar Y1 = pf1.Y()/sp.Y();
    const scalar Y2 = pf2.Y()/sp.Y();

    return BKW<Specie>
    (
        sp,
        Y2*pf2.k_ - Y1*pf1.k_,
        Y2*pf2.kappa_ - Y1*pf1.kappa,
        Y2*pf2.Theta_ - Y1*pf1.Theta_,
        Y2*pf2.a_ - Y1*pf1.a_,
        Y2*pf2.beta_ - Y1*pf1.beta_
    );
}


// ************************************************************************* //
