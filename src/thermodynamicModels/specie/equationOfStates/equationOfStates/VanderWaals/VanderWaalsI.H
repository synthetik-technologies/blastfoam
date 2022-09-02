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

#include "cubicEqn.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::VanderWaals<Specie>::VanderWaals
(
    const Specie& sp,
    const scalar Tc,
    const scalar Pc
)
:
    Specie(sp),
    a_(27.0*sqr(Tc*this->R())/(64.0*Pc)),
    b_(this->R()*Tc/8.0/Pc),
    Tc_(Tc),
    Pc_(Pc)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::VanderWaals<Specie>::VanderWaals
(
    const word& name,
    const VanderWaals<Specie>& vw
)
:
    Specie(name, vw),
    a_(vw.a_),
    b_(vw.b_),
    Tc_(vw.Tc_),
    Pc_(vw.Pc_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::VanderWaals<Specie>>
Foam::VanderWaals<Specie>::clone() const
{
    return autoPtr<VanderWaals<Specie>>
    (
        new VanderWaals<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::VanderWaals<Specie>>
Foam::VanderWaals<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<VanderWaals<Specie>>
    (
        new VanderWaals<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
Foam::scalar Foam::VanderWaals<Specie>::p
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
Foam::scalar Foam::VanderWaals<Specie>::pRhoT
(
    const scalar rho,
    const scalar T
) const
{
    return rho*this->R()*T/(1.0 - b_*rho) - a_*sqr(rho);
}


template<class Specie>
Foam::scalar Foam::VanderWaals<Specie>::rho
(
    const scalar p,
    const scalar T
) const
{
    return cubicEqn
    (
        a_*b_,
        -a_,
        this->R()*T + p*b_,
        p
    ).roots()[0];
}


template<class Specie>
Foam::scalar Foam::VanderWaals<Specie>::Gamma
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
Foam::scalar Foam::VanderWaals<Specie>::cSqr
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
Foam::scalar Foam::VanderWaals<Specie>::dpdv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    const scalar v = 1.0/max(rho, 1e-10);
    return -this->R()*T/sqr(v - b_) + a_*2.0/pow3(v);
}


template<class Specie>
Foam::scalar Foam::VanderWaals<Specie>::dpde
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
Foam::scalar Foam::VanderWaals<Specie>::dpdT
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return this->R()*rho/(1.0 - b_*rho);
}


template<class Specie>
Foam::scalar Foam::VanderWaals<Specie>::E
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return -a_*rho;
}


template<class Specie>
Foam::scalar Foam::VanderWaals<Specie>::Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::VanderWaals<Specie>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::VanderWaals<Specie>::H
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return this->E(rho, e, T) + this->R()*T*b_*rho/(1.0 - b_*rho) - a_*rho;
}


template<class Specie>
Foam::scalar Foam::VanderWaals<Specie>::CpMCv
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar CpCv
) const
{
    return this->R()/(1.0 - b_*rho);
}


template<class Specie>
Foam::scalar Foam::VanderWaals<Specie>::S
(
    const scalar p,
    const scalar T
) const
{
    const scalar rho = this->rho(p, T);
    return
        this->R()
       *(
          - log(p/Pstd)
          + (
              - log(1.0 - b_*rho)/(b_*rho) - 1
              + log(1.0/(1 - b_*rho) - a_*rho/this->R()*T)
            )
        );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::VanderWaals<Specie>::operator+=
(
    const VanderWaals<Specie>& vw
)
{
    scalar X1 = this->Y()/this->W();
    Specie::operator+=(vw);

    if (mag(this->Y()) > small)
    {
        X1 *= this->W()/this->Y();
        const scalar X2 = this->W()*vw.Y()/(vw.W()*this->Y());

        Tc_ = X1*Tc_ + X2*vw.Tc_;
        scalar Vc = 3.0*(X1*b_ + X2*vw.b_);
        Pc_ = 3.0/8.0*RR*Tc_/Vc;
        a_ = 27.0*sqr(Tc_*this->R())/(64.0*Pc_);
        b_ = this->R()*Tc_/8.0/Pc_;
    }
}


template<class Specie>
inline void Foam::VanderWaals<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::VanderWaals<Specie> Foam::operator+
(
    const VanderWaals<Specie>& vw1,
    const VanderWaals<Specie>& vw2
)
{
    Specie sp
    (
        static_cast<const Specie&>(vw1)
      + static_cast<const Specie&>(vw2)
    );

    if (mag(sp.Y()) < small)
    {
        return VanderWaals<Specie>
        (
            sp,
            vw1.Tc_,
            vw1.Pc_
        );
    }
    else
    {
        const scalar X1 = sp.W()*vw1.Y()/(vw1.W()*sp.Y());
        const scalar X2 = sp.W()*vw2.Y()/(vw2.W()*sp.Y());

        const scalar Tc = X1*vw1.Tc_ + X2*vw2.Tc_;
        const scalar Vc = 3.0*(X1*vw1.b_ + X2*vw2.b_);

        return VanderWaals<Specie>
        (
            sp,
            Tc,
            3.0/8.0*sp.R()*Tc/Vc
        );
    }
}


template<class Specie>
inline Foam::VanderWaals<Specie> Foam::operator*
(
    const scalar s,
    const VanderWaals<Specie>& vw
)
{
    return VanderWaals<Specie>
    (
        s*static_cast<const Specie&>(vw),
        vw.Tc_,
        vw.Pc_
    );
}


template<class Specie>
inline Foam::VanderWaals<Specie> Foam::operator==
(
    const VanderWaals<Specie>& vw1,
    const VanderWaals<Specie>& vw2
)
{
    Specie sp
    (
        static_cast<const Specie&>(vw1)
     == static_cast<const Specie&>(vw2)
    );

    const scalar X1 = sp.W()*vw1.Y()/(vw1.W()*sp.Y());
    const scalar X2 = sp.W()*vw2.Y()/(vw2.W()*sp.Y());

    const scalar Tc = X2*vw2.Tc_ - X1*vw1.Tc_;
    const scalar Vc = 3.0*(X2*vw2.b_ - X2*vw1.b_);

    return VanderWaals<Specie>
    (
        sp,
        Tc,
        3.0/8.0*sp.R()*Tc/Vc
    );
}


// ************************************************************************* //
