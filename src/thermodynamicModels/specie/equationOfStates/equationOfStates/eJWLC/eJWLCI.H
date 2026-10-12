/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
inline Foam::scalar Foam::eJWLC<Specie>::V(const scalar rho) const
{
    return rho0_ / max(rho, small);
}


template<class Specie>
inline Foam::scalar Foam::eJWLC<Specie>::ps(const scalar V) const
{
    return
        A_*exp(-R1_*V)
      + B_*exp(-R2_*V)
      + C_*pow(V, -(omega_ + 1.0));
}


template<class Specie>
inline Foam::scalar Foam::eJWLC<Specie>::es(const scalar V) const
{
    // e_s(V) = v0 * [A/R1*exp(-R1*V) + B/R2*exp(-R2*V) + C/omega*V^{-omega}] + E0
    const scalar v0 = 1.0/rho0_;
    return v0*(
        (A_/R1_)*exp(-R1_*V)
      + (B_/R2_)*exp(-R2_*V)
      + (C_/omega_)*pow(V, -omega_)
    ) + E0_;
}


template<class Specie>
inline Foam::scalar Foam::eJWLC<Specie>::Ts(const scalar V) const
{
    // T_s(V) = T_CJ * (V_CJ / V)^omegaT
    const scalar V_CJ = rho0_/max(rho_CJ_, 1e-10);
    return T_CJ_*pow(V_CJ/max(V, 1e-10), omegaT_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::eJWLC<Specie>::eJWLC
(
    const Specie& sp,
    const scalar rho0,
    const scalar omega,
    const scalar A,
    const scalar B,
    const scalar C,
    const scalar R1,
    const scalar R2,
    const scalar T_CJ,
    const scalar rho_CJ,
    const scalar Cv,
    const scalar E0,
    const scalar omegaT
)
:
    Specie(sp),
    rho0_(rho0),
    omega_(omega),
    A_(A),
    B_(B),
    C_(C),
    R1_(R1),
    R2_(R2),
    T_CJ_(T_CJ),
    rho_CJ_(rho_CJ),
    Cv_(Cv),
    E0_(E0),
    omegaT_(omegaT)
{}


template<class Specie>
inline Foam::eJWLC<Specie>::eJWLC
(
    const word& name,
    const eJWLC<Specie>& eos
)
:
    Specie(name, eos),
    rho0_(eos.rho0_),
    omega_(eos.omega_),
    A_(eos.A_),
    B_(eos.B_),
    C_(eos.C_),
    R1_(eos.R1_),
    R2_(eos.R2_),
    T_CJ_(eos.T_CJ_),
    rho_CJ_(eos.rho_CJ_),
    Cv_(eos.Cv_),
    E0_(eos.E0_),
    omegaT_(eos.omegaT_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::eJWLC<Specie>>
Foam::eJWLC<Specie>::clone() const
{
    return autoPtr<eJWLC<Specie>>(new eJWLC<Specie>(*this));
}


template<class Specie>
inline Foam::autoPtr<Foam::eJWLC<Specie>>
Foam::eJWLC<Specie>::New(const dictionary& dict)
{
    return autoPtr<eJWLC<Specie>>(new eJWLC<Specie>(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
Foam::scalar Foam::eJWLC<Specie>::p
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const bool limit
) const
{
    const scalar Vv  = V(rho);
    const scalar pMG = ps(Vv) + omega_*rho*(e - es(Vv));
    return limit ? max(pMG, 0.0) : pMG;
}


template<class Specie>
Foam::scalar Foam::eJWLC<Specie>::pRhoT
(
    const scalar rho,
    const scalar T
) const
{
    // On-isentrope pressure only (T unused — for compatibility)
    return ps(V(rho));
}


template<class Specie>
Foam::scalar Foam::eJWLC<Specie>::rho
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class Specie>
Foam::scalar Foam::eJWLC<Specie>::ePRho
(
    const scalar p,
    const scalar rho
) const
{
    // Invert: p = ps(V) + omega*rho*(e - es(V))  =>  e = es(V) + (p - ps(V)) / (omega*rho)
    const scalar Vv = V(rho);
    return es(Vv) + (p - ps(Vv)) / (omega_*max(rho, 1e-10));
}


template<class Specie>
Foam::scalar Foam::eJWLC<Specie>::Gamma
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar cv
) const
{
    return omega_;
}


template<class Specie>
Foam::scalar Foam::eJWLC<Specie>::cSqr
(
    const scalar p,
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar cv
) const
{
    const scalar rhos = max(rho, small);
    const scalar Vv   = max(V(rhos), 0.1);

    // Exact isentropic sound speed squared:
    // c_s^2 = -(V^2/rho0) * dp_s/dV
    const scalar csSqr =
        (
            A_*R1_*sqr(Vv)*exp(-R1_*Vv)
          + B_*R2_*sqr(Vv)*exp(-R2_*Vv)
          + C_*(omega_ + 1.0)*pow(Vv, -omega_)
        )/rho0_;

    // Thermal (off-isentrope) contribution: omega*(omega + 1)*(e - es(V))
    const scalar cSqrTot = csSqr + omega_*(omega_ + 1.0)*(e - es(Vv));

    return max(cSqrTot, small);
}


template<class Specie>
Foam::scalar Foam::eJWLC<Specie>::dpdv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    // dp/dv = -rho^2 * dp/drho|_e
    const scalar rhos = max(rho, 1e-10);
    const scalar Vv   = V(rhos);

    const scalar dpsdV =
      - A_*R1_*exp(-R1_*Vv)
      - B_*R2_*exp(-R2_*Vv)
      - C_*(omega_ + 1.0)*pow(Vv, -(omega_ + 2.0));

    // dp/dv|_e = dp_s/dv|_e + omega*rho * (-des/dv|_e) - omega*(e - es)*rho^2
    // dp_s/dv = dp_s/dV * dV/dv = dpsdV * rho0 (since V = rho0*v)
    // des/dv  = des/dV * rho0 = -ps(V)*rho0/rho0 = -ps(V) ... wait
    // v = 1/rho, V = rho0/rho = rho0*v, dV/dv = rho0
    const scalar dpsdv_e = dpsdV * rho0_;

    // des/dv = des/dV * dV/dv = (-ps(V)/rho0) * rho0 = -ps(V)
    const scalar desdv_e = -ps(Vv);

    // MG off-isentrope: omega*rho*(e - es(V)), d/dv at const e:
    // d/dv [omega*rho*(e-es)] = omega*(e-es)*d(rho)/dv + omega*rho*(-des/dv)
    //                         = omega*(e-es)*(-rho^2) + omega*rho*ps(V)
    return dpsdv_e
         + omega_*(e - es(Vv))*(-sqr(rhos))
         + omega_*rhos*(-desdv_e);
}


template<class Specie>
Foam::scalar Foam::eJWLC<Specie>::dpde
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return omega_*rho;
}


template<class Specie>
Foam::scalar Foam::eJWLC<Specie>::dpdT
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    // Zero for energy-based EOS — temperature enters only through the
    // caloric model, not through the EOS pressure directly.
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::eJWLC<Specie>::E
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    const scalar Vv = V(max(rho, 1e-10));
    return es(Vv) - Cv_*Ts(Vv);
}


template<class Specie>
Foam::scalar Foam::eJWLC<Specie>::H
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    // Enthalpy correction: H = E + p*v = E + p/rho
    const scalar rhos = max(rho, 1e-10);
    const scalar Vv   = V(rhos);
    return es(Vv) - Cv_*Ts(Vv) + ps(Vv)/rhos;
}


template<class Specie>
Foam::scalar Foam::eJWLC<Specie>::Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    // The EOS layer contributes zero to Cv; the caloric model provides it.
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::eJWLC<Specie>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::eJWLC<Specie>::CpMCv
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar CpCv
) const
{
    // Cp - Cv = Gamma * Cv   (from Mie-Grüneisen relation, eBased)
    return this->eBased_ ? omega_*CpCv : CpCv*omega_/(omega_ + 1.0);
}


template<class Specie>
Foam::scalar Foam::eJWLC<Specie>::Sp
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::eJWLC<Specie>::Sv
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
inline void Foam::eJWLC<Specie>::operator+=
(
    const eJWLC<Specie>& eos
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(eos);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = eos.Y()/this->Y();

        rho0_   = Y1*rho0_   + Y2*eos.rho0_;
        omega_  = Y1*omega_  + Y2*eos.omega_;
        A_      = Y1*A_      + Y2*eos.A_;
        B_      = Y1*B_      + Y2*eos.B_;
        C_      = Y1*C_      + Y2*eos.C_;
        R1_     = Y1*R1_     + Y2*eos.R1_;
        R2_     = Y1*R2_     + Y2*eos.R2_;
        T_CJ_   = Y1*T_CJ_  + Y2*eos.T_CJ_;
        rho_CJ_ = Y1*rho_CJ_+ Y2*eos.rho_CJ_;
        Cv_     = Y1*Cv_     + Y2*eos.Cv_;
        E0_     = Y1*E0_     + Y2*eos.E0_;
        omegaT_ = Y1*omegaT_ + Y2*eos.omegaT_;
    }
}


template<class Specie>
inline void Foam::eJWLC<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::eJWLC<Specie> Foam::operator+
(
    const eJWLC<Specie>& eos1,
    const eJWLC<Specie>& eos2
)
{
    Specie sp
    (
        static_cast<const Specie&>(eos1)
      + static_cast<const Specie&>(eos2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return eJWLC<Specie>
        (
            sp,
            eos1.rho0_,
            eos1.omega_,
            eos1.A_,
            eos1.B_,
            eos1.C_,
            eos1.R1_,
            eos1.R2_,
            eos1.T_CJ_,
            eos1.rho_CJ_,
            eos1.Cv_,
            eos1.E0_,
            eos1.omegaT_
        );
    }
    else
    {
        const scalar Y1 = eos1.Y()/sp.Y();
        const scalar Y2 = eos2.Y()/sp.Y();

        return eJWLC<Specie>
        (
            sp,
            Y1*eos1.rho0_    + Y2*eos2.rho0_,
            Y1*eos1.omega_   + Y2*eos2.omega_,
            Y1*eos1.A_       + Y2*eos2.A_,
            Y1*eos1.B_       + Y2*eos2.B_,
            Y1*eos1.C_       + Y2*eos2.C_,
            Y1*eos1.R1_      + Y2*eos2.R1_,
            Y1*eos1.R2_      + Y2*eos2.R2_,
            Y1*eos1.T_CJ_    + Y2*eos2.T_CJ_,
            Y1*eos1.rho_CJ_  + Y2*eos2.rho_CJ_,
            Y1*eos1.Cv_      + Y2*eos2.Cv_,
            Y1*eos1.E0_      + Y2*eos2.E0_,
            Y1*eos1.omegaT_  + Y2*eos2.omegaT_
        );
    }
}


template<class Specie>
inline Foam::eJWLC<Specie> Foam::operator*
(
    const scalar s,
    const eJWLC<Specie>& eos
)
{
    return eJWLC<Specie>
    (
        s*static_cast<const Specie&>(eos),
        eos.rho0_,
        eos.omega_,
        eos.A_,
        eos.B_,
        eos.C_,
        eos.R1_,
        eos.R2_,
        eos.T_CJ_,
        eos.rho_CJ_,
        eos.Cv_,
        eos.E0_,
        eos.omegaT_
    );
}


template<class Specie>
inline Foam::eJWLC<Specie> Foam::operator==
(
    const eJWLC<Specie>& eos1,
    const eJWLC<Specie>& eos2
)
{
    Specie sp
    (
        static_cast<const Specie&>(eos1)
     == static_cast<const Specie&>(eos2)
    );

    const scalar Y1 = eos1.Y()/sp.Y();
    const scalar Y2 = eos2.Y()/sp.Y();

    return eJWLC<Specie>
    (
        sp,
        Y2*eos2.rho0_    - Y1*eos1.rho0_,
        Y2*eos2.omega_   - Y1*eos1.omega_,
        Y2*eos2.A_       - Y1*eos1.A_,
        Y2*eos2.B_       - Y1*eos1.B_,
        Y2*eos2.C_       - Y1*eos1.C_,
        Y2*eos2.R1_      - Y1*eos1.R1_,
        Y2*eos2.R2_      - Y1*eos1.R2_,
        Y2*eos2.T_CJ_    - Y1*eos1.T_CJ_,
        Y2*eos2.rho_CJ_  - Y1*eos1.rho_CJ_,
        Y2*eos2.Cv_      - Y1*eos1.Cv_,
        Y2*eos2.E0_      - Y1*eos1.E0_,
        Y2*eos2.omegaT_  - Y1*eos1.omegaT_
    );
}


// ************************************************************************* //
