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
inline Foam::Tillotson<Specie>::Tillotson
(
    const Specie& sp,
    const scalar rho0,
    const scalar e0,
    const scalar a,
    const scalar b,
    const scalar A,
    const scalar B,
    const scalar pCav,
    const scalar rhoCav,
    const scalar k
)
:
    Specie(sp),
    rho0_(rho0),
    e0_(e0),
    a_(a),
    b_(b),
    A_(A),
    B_(B),
    pCav_(pCav),
    rhoCav_(rhoCav),
    k_(k)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::Tillotson<Specie>::Tillotson
(
    const word& name,
    const Tillotson<Specie>& pf
)
:
    Specie(name, pf),
    rho0_(pf.rho0_),
    e0_(pf.e0_),
    a_(pf.a_),
    b_(pf.b_),
    A_(pf.A_),
    B_(pf.B_),
    pCav_(pf.pCav_),
    rhoCav_(pf.rhoCav_),
    k_(pf.k_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::Tillotson<Specie>>
Foam::Tillotson<Specie>::clone() const
{
    return autoPtr<Tillotson<Specie>>
    (
        new Tillotson<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::Tillotson<Specie>>
Foam::Tillotson<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<Tillotson<Specie>>
    (
        new Tillotson<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::p
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const bool limit
) const
{
    scalar eta(rho/rho0_);
    scalar mu(eta - 1.0);
    scalar res =
        (a_ + b_/(e/(e0_*sqr(eta)) + 1.0))*rho*e
      + A_*mu
      + B_*sqr(mu);

      return limit ? max(res, pCav_) : res;
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::Pi
(
    const scalar rho,
    const scalar e
) const
{
    if (rho < rhoCav_)
    {
        return 0.0;
    }

    scalar mu(rho/rho0_ - 1.0);
    return
      - A_*mu
      - B_*sqr(mu)
      - b_*sqr(e*rho*rho0_)*e0_*rho
       /sqr(e*sqr(rho0_) + e0_*sqr(rho));
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::delta
(
    const scalar p,
    const scalar rho,
    const scalar e
) const
{
    if (rho < rhoCav_)
    {
        return 0.0;
    }

    scalar x(e*sqr(rho0_) + e0_*sqr(rho));
    scalar G(GammaMG(rho, e));
    scalar dGammadRho =
        4.0*b_*e*sqr(e0_*rho0_)*pow3(rho)/pow3(x);

    scalar dPidRho =
      - A_/rho0_
      - B_*2.0*(rho - rho0_)/sqr(rho0_)
      + b_*e*e0_*sqr(rho*rho0_)/pow3(x)*(e0_*sqr(rho) - 3.0*e*sqr(rho0_));
    return
      - (p + Pi(rho, e))*dGammadRho/sqr(G - 1.0)
      + dPidRho/(G - 1.0);
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::cSqr
(
    const scalar p,
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar Cv
) const
{
//     if (rho < rhoCav_)
//     {
//         return sqrt(k_*pCav_/max(rho, 1e-10));
//     }

    scalar gamma(GammaMG(rho, e));
    scalar h((gamma*p + Pi(rho, e))/((gamma - 1.0)*max(rho, 1e-10)));
    scalar cLSqr((h - delta(p, rho, e))*(gamma - 1.0));

//     if (rho < rho0_)
//     {
//         scalar x((rho - rho0_)/(rhoCav_ - rho0_));
//         scalar cCavSqr(k_*pCav_/max(rho, 1e-10));
//
//         return
//             1.0
//            /(
//                 max(rho, 1e-10)
//                *(
//                     x/(max(rho, 1e-10)*cCavSqr)
//                   + (1.0 - x)/(rho0_*cLSqr)
//                 )
//             );
//     }
    return cLSqr;
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::GammaMG
(
    const scalar rho,
    const scalar e
) const
{
//     if (rho < rhoCav_)
//     {
//         return k_;
//     }
    scalar x(e*sqr(rho0_) + e0_*sqr(rho));
    scalar w =
        (a_*x + b_*e0_*sqr(rho))/x
      - b_*e*e0_*sqr(rho*rho0_)/sqr(x);
//     else if (rho < rho0_)
//     {
//         scalar x((rho - rho0_)/(rhoCav_ - rho0_));
//
//         return
//             1.0
//            /(
//                 max(rho, 1e-10)
//                *(
//                     x/(rhoCav_*k_)
//                   + (1.0 - x)/(rho0_*(omega_ + 1.0))
//                 )
//             );
//     }

    return w + 1.0;
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::dpdv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    scalar x(e*sqr(rho0_) + e0_*sqr(rho));
    return
       -sqr(rho/rho0_*x)
       *(
            A_*rho0_*sqr(x)
          + 2.0*B_*sqr(x)*(rho - rho0_)
          + 2.0*b_*sqr(e*rho*sqr(rho0_))*e0_
          + e*sqr(rho0_)*x*(a_*x + b_*e0_*sqr(rho))
        );
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::dpde
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    scalar x(e*sqr(rho0_) + e0_*sqr(rho));
    return
        rho/sqr(x)
       *(
          - b_*e*e0_*sqr(rho*rho0_)
          + x*(a_*x + b_*e0_*sqr(rho))
        );
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::dpdT
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    NotImplemented;
    return this->R()*rho;
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::E
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return EcTable_.lookup(rho);
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    // Check
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    // Check
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::H
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::S
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
inline void Foam::Tillotson<Specie>::operator+=
(
    const Tillotson<Specie>& pf
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(pf);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = pf.Y()/this->Y();

        rho0_ = Y1*rho0_ + Y2*pf.rho0_;
        e0_ = Y1*e0_ + Y2*pf.e0_;
        a_ = Y1*a_ + Y2*pf.a_;
        b_ = Y1*b_ + Y2*pf.b_;
        A_ = Y1*A_ + Y2*pf.A_;
        B_ = Y1*B_ + Y2*pf.B_;
        pCav_ = Y1*pCav_ + Y2*pf.pCav_;
        rhoCav_ = Y1*rhoCav_ + Y2*pf.rhoCav_;
        k_ = Y1*k_ + Y2*pf.k_;
    }
}


template<class Specie>
inline void Foam::Tillotson<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::Tillotson<Specie> Foam::operator+
(
    const Tillotson<Specie>& pf1,
    const Tillotson<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
      + static_cast<const Specie&>(pf2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return Tillotson<Specie>
        (
            sp,
            pf1.rho0_,
            pf1.e0_,
            pf1.a_,
            pf1.b_,
            pf1.A_,
            pf1.B_,
            pf1.pCav_,
            pf1.rhoCav_,
            pf1.k_
        );
    }
    else
    {
        const scalar Y1 = pf1.Y()/sp.Y();
        const scalar Y2 = pf2.Y()/sp.Y();

        return Tillotson<Specie>
        (
            sp,
            Y1*pf1.rho0_ + Y2*pf2.rho0_,
            Y1*pf1.e0_ + Y2*pf2.e0_,
            Y1*pf1.a_ + Y2*pf2.a_,
            Y1*pf1.b_ + Y2*pf2.b_,
            Y1*pf1.A_ + Y2*pf2.A_,
            Y1*pf1.B_ + Y2*pf2.B_,
            Y1*pf1.pCav_ + Y2*pf2.pCav_,
            Y1*pf1.rhoCav_ + Y2*pf2.rhoCav_,
            Y1*pf1.k_ + Y2*pf2.k_
        );
    }
}


template<class Specie>
inline Foam::Tillotson<Specie> Foam::operator*
(
    const scalar s,
    const Tillotson<Specie>& pf
)
{
    return Tillotson<Specie>
    (
        s*static_cast<const Specie&>(pf),
        pf.rho0_,
        pf.e0_,
        pf.a_,
        pf.b_,
        pf.A_,
        pf.B_,
        pf.pCav_,
        pf.rhoCav_,
        pf.k_
    );
}


template<class Specie>
inline Foam::Tillotson<Specie> Foam::operator==
(
    const Tillotson<Specie>& pf1,
    const Tillotson<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
     == static_cast<const Specie&>(pf2)
    );

    const scalar Y1 = pf1.Y()/sp.Y();
    const scalar Y2 = pf2.Y()/sp.Y();

    return Tillotson<Specie>
    (
        sp,
        Y2*pf2.rho0_    - Y1*pf1.rho0_,
        Y2*pf2.e0_    - Y1*pf1.e0_,
        Y2*pf2.a_  - Y1*pf1.a_,
        Y2*pf2.b_ - Y1*pf1.b_,
        Y2*pf2.A_  - Y1*pf1.A_,
        Y2*pf2.B_ - Y1*pf1.B_,
        Y2*pf2.pCav_  - Y1*pf1.pCav_,
        Y2*pf2.rhoCav_  - Y1*pf1.rhoCav_,
        Y2*pf2.k_    - Y1*pf1.k_
    );
}


// ************************************************************************* //
