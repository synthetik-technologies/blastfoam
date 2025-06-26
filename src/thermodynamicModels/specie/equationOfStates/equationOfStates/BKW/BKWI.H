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
    scalar dxdrho = -a_*kappa_*k_/(this->W()*pow(T + Theta_, a_ + 1.0));
    return
      - sqr(rho)*this->R()*T
       *(
            1.0
          + exp(beta_*x)*(x + rho*dxdrho*(1.0 + x*beta_))
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
    return rho*this->R()*(1.0 + exp(beta_*x)*(x + T*dxdT*(1.0 + x*beta_)));
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
    scalar dxdT = -a_*kappa_*k_*rho/(this->W()*pow(T + Theta_, a_ + 1.0));
    return
        this->R()*T*a_*exp(beta_*x)/(T + Theta_)
       *(
            x*(2.0 - T/(T + Theta_))
          + T*dxdT*(1.0 + beta_*x)
        );
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
Foam::scalar Foam::BKW<Specie>::Sp
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::BKW<Specie>::Sv
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
    scalar X1 = this->Y()/this->W();
    Specie::operator+=(pf);

    if (mag(this->Y()) > SMALL)
    {
        X1 *= this->W()/this->Y();
        const scalar X2 = this->W()*pf.Y()/(pf.W()*this->Y());

        k_ = X1*k_ + X2*pf.k_;
        kappa_ = X1*kappa_ + X2*pf.kappa_;
        Theta_ = X1*Theta_ + X2*pf.Theta_;
        a_ = X1*a_ + X2*pf.a_;
        beta_ = X1*beta_ + X2*pf.beta_;
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
        const scalar X1 = sp.W()*pf1.Y()/(pf1.W()*sp.Y());
        const scalar X2 = sp.W()*pf2.Y()/(pf2.W()*sp.Y());

        return BKW<Specie>
        (
            sp,
            X1*pf1.b_ + X2*pf2.b_
        );

        return BKW<Specie>
        (
            sp,
            X1*pf1.k_ + X2*pf2.k_,
            X1*pf1.kappa_ + X2*pf2.kappa,
            X1*pf1.Theta_ + X2*pf2.Theta_,
            X1*pf1.a_ + X2*pf2.a_,
            X1*pf1.beta_ + X2*pf2.beta_
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
