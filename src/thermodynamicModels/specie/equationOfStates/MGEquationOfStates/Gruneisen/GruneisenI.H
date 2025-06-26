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
inline Foam::Gruneisen<Specie>::Gruneisen
(
    const Specie& sp,
    const scalar rho0,
    const scalar C,
    const scalar S1,
    const scalar S2,
    const scalar S3,
    const scalar Gamma,
    const scalar a
)
:
    Specie(sp),
    rho0_(rho0),
    C_(C),
    S1_(S1),
    S2_(S2),
    S3_(S3),
    Gamma_(Gamma),
    a_(a)
{}


template<class Specie>
Foam::scalar Foam::Gruneisen<Specie>::Pi
(
    const scalar rho,
    const scalar e
) const
{
    scalar mu(rho/rho0_ - 1.0);
    if (mu <= 0)
    {
        return -rho0_*sqr(C_)*mu;
    }
    return
      - rho0_*sqr(C_)*mu
       *(
            1.0
          + (1.0 - 0.5*Gamma_)*mu
          - 0.5*a_*sqr(mu)
        )
       /(
            1.0
          - (S1_ - 1.0)*mu
          - S2_*sqr(mu)/(mu + 1.0)
          - S3_*pow3(mu)/sqr(mu + 1.0)
        );
}


template<class Specie>
Foam::scalar Foam::Gruneisen<Specie>::dPidRho
(
    const scalar rho,
    const scalar e
) const
{
    scalar mu(rho/rho0_ - 1.0);
    if (rho <= rho0_)
    {
        return -2.0*sqr(C_);
    }

    scalar dmudRho(1.0/rho0_);
    scalar num =
        mu
      + (1.0 - 0.5*Gamma_)*sqr(mu)
      - 0.5*a_*pow3(mu);
    scalar dnum =
        (
            1.0
          + (2.0 - Gamma_)*mu
          - 1.5*a_*sqr(mu)
        )/rho0_;
    scalar denom =
        1.0
      - (S1_ - 1.0)*mu
      - S2_*sqr(mu)/(mu + 1.0)
      - S3_*pow3(mu)/sqr(mu + 1.0);
    scalar ddenom =
        (
          - (S1_ - 1.0)
          - S2_*mu/(mu + 1.0)*(2.0 - mu/(mu + 1))
          - S3_*sqr(mu/(mu + 1.0))*(3.0 - 2.0*mu/sqr(mu + 1))
        )*dmudRho;
    return
        -rho0_*sqr(C_)*(dnum/sqr(denom) - 2.0*num/pow3(denom)*ddenom);
}


template<class Specie>
Foam::scalar Foam::Gruneisen<Specie>::delta
(
    const scalar p,
    const scalar rho,
    const scalar e
) const
{
    return
      - (p + Pi(rho, e))*dGammaMGdRho(rho, e)
       /sqr(GammaMG(rho, e) - 1.0)
     + dPidRho(rho, e)/(GammaMG(rho, e) - 1.0);
}


template<class Specie>
Foam::scalar Foam::Gruneisen<Specie>::GammaMG
(
    const scalar rho,
    const scalar e
) const
{
    return Gamma_ + a_*(rho/rho0_ - 1.0) + 1.0;
}

template<class Specie>
Foam::scalar Foam::Gruneisen<Specie>::dGammaMGdRho
(
    const scalar rho,
    const scalar e
) const
{
    return a_/rho0_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::Gruneisen<Specie>::Gruneisen
(
    const word& name,
    const Gruneisen<Specie>& pf
)
:
    Specie(name, pf),
    rho0_(pf.rho0_),
    C_(pf.C_),
    S1_(pf.S1_),
    S2_(pf.S2_),
    S3_(pf.S3_),
    Gamma_(pf.Gamma_),
    a_(pf.a_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::Gruneisen<Specie>>
Foam::Gruneisen<Specie>::clone() const
{
    return autoPtr<Gruneisen<Specie>>
    (
        new Gruneisen<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::Gruneisen<Specie>>
Foam::Gruneisen<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<Gruneisen<Specie>>
    (
        new Gruneisen<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
Foam::scalar Foam::Gruneisen<Specie>::p
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const bool limit
) const
{
    scalar pMG = (GammaMG(rho, e) - 1.0)*rho*e - Pi(rho, e);
    return limit ? max(pMG, 0.0) : pMG;
}


template<class Specie>
Foam::scalar Foam::Gruneisen<Specie>::rho
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class Specie>
Foam::scalar Foam::Gruneisen<Specie>::Gamma
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar Cv
) const
{
    return GammaMG(rho, e);
}


template<class Specie>
Foam::scalar Foam::Gruneisen<Specie>::cSqr
(
    const scalar p,
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar Cv
) const
{
    scalar G = GammaMG(rho, e);
    scalar h((G*p + Pi(rho, e))/((G - 1.0)*max(rho, 1e-10)));
    return (h - delta(p, rho, e))*(G - 1.0);
}


template<class Specie>
Foam::scalar Foam::Gruneisen<Specie>::dpdv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return
      - sqr(rho)
       *(
           GammaMG(rho, e)*e
         + dGammaMGdRho(rho, e)*rho*e
         - dPidRho(rho, e)
        );
}


template<class Specie>
Foam::scalar Foam::Gruneisen<Specie>::dpde
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return GammaMG(rho, e)*rho;
}


template<class Specie>
Foam::scalar Foam::Gruneisen<Specie>::dpdT
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
Foam::scalar Foam::Gruneisen<Specie>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::Gruneisen<Specie>::Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::Gruneisen<Specie>::CpMCv
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar CpCv
) const
{
    scalar G = GammaMG(rho, e);
    return this->eBased_ ? (G - 1.0)*CpCv : CpCv*(G - 1.0)/G;
}


template<class Specie>
Foam::scalar Foam::Gruneisen<Specie>::E
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::Gruneisen<Specie>::H
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::Gruneisen<Specie>::Sp
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::Gruneisen<Specie>::Sv
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
inline void Foam::Gruneisen<Specie>::operator+=
(
    const Gruneisen<Specie>& pf
)
{
    NotImplemented;
//     scalar Y1 = this->Y();
//     Specie::operator+=(pf);
//
//     if (mag(this->Y()) > SMALL)
//     {
//         Y1 /= this->Y();
//         const scalar Y2 = pf.Y()/this->Y();
//
//         gamma_ = Y1*gamma_ + Y2*pf.gamma_;
//     }
}


template<class Specie>
inline void Foam::Gruneisen<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::Gruneisen<Specie> Foam::operator+
(
    const Gruneisen<Specie>& pf1,
    const Gruneisen<Specie>& pf2
)
{
    NotImplemented;
    return pf1;
//     Specie sp
//     (
//         static_cast<const Specie&>(pf1)
//       + static_cast<const Specie&>(pf2)
//     );
//
//     if (mag(sp.Y()) < SMALL)
//     {
//         return Gruneisen<Specie>
//         (
//             sp,
//             pf1.gamma_
//         );
//     }
//     else
//     {
//         const scalar Y1 = pf1.Y()/sp.Y();
//         const scalar Y2 = pf2.Y()/sp.Y();
//
//         return Gruneisen<Specie>
//         (
//             sp,
//             Y1*pf1.gamma_ + Y2*pf2.gamma_
//         );
//     }
}


template<class Specie>
inline Foam::Gruneisen<Specie> Foam::operator*
(
    const scalar s,
    const Gruneisen<Specie>& pf
)
{

    NotImplemented;
    return pf;
//     return Gruneisen<Specie>
//     (
//         s*static_cast<const Specie&>(pf),
//         pf.gamma_
//     );
}


template<class Specie>
inline Foam::Gruneisen<Specie> Foam::operator==
(
    const Gruneisen<Specie>& pf1,
    const Gruneisen<Specie>& pf2
)
{
    NotImplemented;
    return pf1;
//     Specie sp
//     (
//         static_cast<const Specie&>(pf1)
//      == static_cast<const Specie&>(pf2)
//     );
//
//     const scalar Y1 = pf1.Y()/sp.Y();
//     const scalar Y2 = pf2.Y()/sp.Y();
//
//     return Gruneisen<Specie>
//     (
//         sp,
//         Y2*pf2.gamma_ - Y1*pf1.gamma_
//     );
}


// ************************************************************************* //
