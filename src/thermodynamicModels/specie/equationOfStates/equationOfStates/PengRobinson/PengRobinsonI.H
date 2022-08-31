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

#include "mathematicalConstants.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::PengRobinson<Specie>::PengRobinson
(
    const Specie& sp,
    const scalar Tc,
    const scalar Vc,
    const scalar Pc,
    const scalar omega,
    const scalar Zc
)
:
    Specie(sp),
    Tc_(Tc),
    Vc_(Vc),
    Pc_(Pc),
    omega_(omega),
    Zc_(Zc),
    a_(0.45724*sqr(RR*Tc_)/Pc_),
    b_(0.07780*RR*Tc_/Pc_),
    kappa_(0.37464 + 1.54226*omega_ - 0.26992*sqr(omega_))
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::PengRobinson<Specie>::PengRobinson
(
    const word& name,
    const PengRobinson<Specie>& pr
)
:
    Specie(name, pr),
    Tc_(pr.Tc_),
    Vc_(pr.Vc_),
    Pc_(pr.Pc_),
    omega_(pr.omega_),
    Zc_(pr.Zc_),
    a_(pr.a_),
    b_(pr.b_),
    kappa_(pr.kappa_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::PengRobinson<Specie>>
Foam::PengRobinson<Specie>::clone() const
{
    return autoPtr<PengRobinson<Specie>>
    (
        new PengRobinson<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::PengRobinson<Specie>>
Foam::PengRobinson<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<PengRobinson<Specie>>
    (
        new PengRobinson<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
Foam::scalar Foam::PengRobinson<Specie>::p
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
Foam::scalar Foam::PengRobinson<Specie>::pRhoT
(
    const scalar rho,
    const scalar T
) const
{
    const scalar Vm = this->W()/max(rho, 1e-10);
    const scalar alpha = sqr(1.0 + kappa_*(1.0 - sqrt(T/Tc_)));

    return RR*T/(Vm - b_) - a_*alpha/(sqr(Vm) + 2.0*b_*Vm - sqr(b_));
}


template<class Specie>
inline Foam::scalar Foam::PengRobinson<Specie>::Z
(
    scalar p,
    scalar T
) const
{
    const scalar Tr = T/Tc_;
    const scalar alpha = sqr(1 + kappa_*(1 - sqrt(Tr)));

    const scalar A = a_*alpha*p/sqr(RR*T);
    const scalar B = b_*p/(RR*T);

    const scalar a2 = B - 1;
    const scalar a1 = A - 2*B - 3*sqr(B);
    const scalar a0 = -A*B + sqr(B) + pow3(B);

    const scalar Q = (3*a1 - a2*a2)/9.0;
    const scalar Rl = (9*a2*a1 - 27*a0 - 2*a2*a2*a2)/54.0;

    const scalar Q3 = Q*Q*Q;
    const scalar D = Q3 + Rl*Rl;

    scalar root = -1;

    if (D <= 0)
    {
        const scalar th = ::acos(Rl/sqrt(-Q3));
        const scalar qm = 2*sqrt(-Q);
        const scalar r1 = qm*cos(th/3.0) - a2/3.0;
        const scalar r2 =
            qm*cos((th + 2*constant::mathematical::pi)/3.0) - a2/3.0;
        const scalar r3 =
            qm*cos((th + 4*constant::mathematical::pi)/3.0) - a2/3.0;

        root = max(r1, max(r2, r3));
    }
    else
    {
        // One root is real
        const scalar D05 = sqrt(D);
        const scalar S = pow(Rl + D05, 1.0/3.0);
        scalar Tl = 0;
        if (D05 > Rl)
        {
            Tl = -pow(mag(Rl - D05), 1.0/3.0);
        }
        else
        {
            Tl = pow(Rl - D05, 1.0/3.0);
        }

        root = S + Tl - a2/3.0;
    }

    return root;
}


template<class Specie>
Foam::scalar Foam::PengRobinson<Specie>::rho
(
    const scalar p,
    const scalar T
) const
{
    return p/(this->Z(p, T)*this->R()*T);
}


template<class Specie>
Foam::scalar Foam::PengRobinson<Specie>::Gamma
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
Foam::scalar Foam::PengRobinson<Specie>::cSqr
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
Foam::scalar Foam::PengRobinson<Specie>::dpdv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    const scalar W = this->W();
    const scalar Vm = W/max(rho, 1e-10);
    const scalar Tr = T/Tc_;
    const scalar alpha = sqr(1 + kappa_*(1.0 - sqrt(Tr)));

    return
      - RR*T*W/sqr(Vm - b_)
      + a_*alpha
       *2.0*W*(Vm + b_)
       /sqr(sqr(Vm) + 2.0*Vm*b_ - sqr(b_));
}


template<class Specie>
Foam::scalar Foam::PengRobinson<Specie>::dpde
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
Foam::scalar Foam::PengRobinson<Specie>::dpdT
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    const scalar Vm = this->W()/max(rho, 1e-10);
    const scalar Tr = T/Tc_;
    const scalar dalphadT = (1.0 + kappa_*(1.0 - sqrt(Tr)))*kappa_/sqrt(Tc_*Tr);

    return RR/(Vm - b_) + a_*dalphadT/(sqr(Vm) + 2.0*b_*Vm - sqr(b_));
}


template<class Specie>
Foam::scalar Foam::PengRobinson<Specie>::E
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    // const scalar Vm = this->W()/max(rho, 1e-10);
    // const scalar sqrtTr = sqrt(T/Tc_);
    // const scalar sqrt2 = sqrt(2.0);
    // return
    //   - a_*(1.0 + kappa_*(1.0 - sqrtTr))
    //    *(
    //         kappa_/Tc_/sqrtTr
    //       - (1.0 + kappa_*(1.0 - sqrtTr))
    //     )/(2.0*sqrt(2.0)*b_*this->W())
    //    *log((Vm + b_*(sqrt2 + 1.0))/(Vm - b_*(sqrt2 - 1.0))); // Double check

    const scalar p = this->pRhoT(rho, T);
    const scalar Pr = p/Pc_;
    const scalar Tr = T/Tc_;
    const scalar B = 0.07780*Pr/Tr;
    const scalar alpha = sqr(1 + kappa_*(1 - sqrt(Tr)));

    const scalar Z = this->Z(p, T);

    return
        this->R()
       *Tc_
       *(
         - 2.078*(1 + kappa_)*sqrt(alpha)
          *log((Z + 2.414*B)/(Z - 0.414*B))
        );
}


template<class Specie>
Foam::scalar Foam::PengRobinson<Specie>::Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    const scalar sqrtTr = sqrt(T/Tc_);
    const scalar sqrt2 = sqrt(2.0);
    return
        a_*kappa_*(1.0 + kappa_)/(2.0*sqrt(Tc_)*pow3(sqrtTr))
       *log((1.0 + (1.0 - sqrt2)*b_*rho)/(1.0 + (1.0 - sqrt2)*b_*rho)) // Double check
       /(2.0*sqrt2*b_*this->W());

    // const scalar p = this->pRhoT(rho, T);
    // const scalar B = b_*p/(RR*T);
    // const scalar Z = this->Z(p, T);

    // const scalar app = kappa_*a_*(1 + kappa_)/(2*sqrt(pow3(T)*Tc_));

    // const scalar root2 = sqrt(2.0);

    // return
    // (
    //     app*(T/(2*root2*b_))*log((Z + (root2 + 1)*B)/(Z - (root2 - 1)*B))
    //   - RR
    // )/this->W();
}


template<class Specie>
Foam::scalar Foam::PengRobinson<Specie>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Cv(rho, e, T) - T*sqr(dpdT(rho, e, T))/dpdv(rho, e, T) - this->R();
}


template<class Specie>
Foam::scalar Foam::PengRobinson<Specie>::H
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    // int(v*dpdv + T*dpdT, v)
    const scalar p = this->pRhoT(rho, T);
    const scalar Pr = p/Pc_;
    const scalar Tr = T/Tc_;
    const scalar B = 0.07780*Pr/Tr;
    const scalar alpha = sqr(1 + kappa_*(1 - sqrt(Tr)));

    const scalar Z = this->Z(p, T);

    return
        this->R()
       *Tc_
       *(
           Tr*(Z - 1)
         - 2.078*(1 + kappa_)*sqrt(alpha)
          *log((Z + 2.414*B)/(Z - 0.414*B))
        );
}


template<class Specie>
Foam::scalar Foam::PengRobinson<Specie>::CpMCv
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar CpCv
) const
{
    return dpdT(rho, e, T)/max(rho, 1e-10);

    // const scalar p = this->pRhoT(rho, T);
    // const scalar Tr = T/Tc_;
    // const scalar alpha = sqr(1 + kappa_*(1 - sqrt(Tr)));

    // const scalar A = alpha*a_*p/sqr(RR*T);
    // const scalar B = b_*p/(RR*T);

    // const scalar Z = this->Z(p, T);

    // const scalar ap = kappa_*a_*(kappa_/Tc_ - (1 + kappa_)/sqrt(T*Tc_));
    // const scalar M = (sqr(Z) + 2*B*Z - sqr(B))/(Z - B);
    // const scalar N = ap*B/(b_*RR);

    // return this->R()*sqr(M - N)/(sqr(M) - 2*A*(Z + B));
}


template<class Specie>
Foam::scalar Foam::PengRobinson<Specie>::S
(
    const scalar p,
    const scalar T
) const
{
    const scalar Pr = p/Pc_;
    const scalar Tr = T/Tc_;
    const scalar B = 0.07780*Pr/Tr;

    const scalar Z = this->Z(p, T);

    return
        this->R()
       *(
          - log(p/Pstd)
          + (
                log(Z - B)
              - 2.078*kappa_*((1 + kappa_)/sqrt(Tr) - kappa_)
               *log((Z + 2.414*B)/(Z - 0.414*B))
            )
        );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::PengRobinson<Specie>::operator+=
(
    const PengRobinson<Specie>& pr
)
{
    scalar X1 = this->Y()/this->W();
    Specie::operator+=(pr);

    if (mag(this->Y()) > small)
    {
        X1 *= this->W()/this->Y();
        const scalar X2 = this->W()*pr.Y()/(pr.W()*this->Y());

        Tc_ = X1*Tc_ + X2*pr.Tc_;
        Vc_ = X1*Vc_ + X2*pr.Vc_;
        Pc_ = RR*Zc_*Tc_/Vc_;
        omega_ = X1*omega_ + X2*pr.omega_;
        Zc_ = X1*Zc_ + X2*pr.Zc_;
        a_ = 0.45724*sqr(RR*Tc_)/Pc_;
        b_ = 0.07780*RR*Tc_/Pc_;
        kappa_ = 0.37464 + 1.54226*omega_ - 0.26992*sqr(omega_);
    }
}


template<class Specie>
inline void Foam::PengRobinson<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::PengRobinson<Specie> Foam::operator+
(
    const PengRobinson<Specie>& pr1,
    const PengRobinson<Specie>& pr2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pr1)
      + static_cast<const Specie&>(pr2)
    );

    if (mag(sp.Y()) < small)
    {
        return PengRobinson<Specie>
        (
            sp,
            pr1.Tc_,
            pr1.Vc_,
            pr1.Pc_,
            pr1.omega_,
            pr1.Zc_
        );
    }
    else
    {
        const scalar X1 = sp.W()*pr1.Y()/(pr1.W()*sp.Y());
        const scalar X2 = sp.W()*pr2.Y()/(pr2.W()*sp.Y());

        const scalar Tc = X1*pr1.Tc_ + X2*pr2.Tc_;
        const scalar Vc = X1*pr1.Vc_ + X2*pr2.Vc_;
        const scalar Zc = X1*pr1.Zc_ + X2*pr2.Zc_;

        return PengRobinson<Specie>
        (
            sp,
            Tc,
            Vc,
            RR*Zc*Tc/Vc,
            X1*pr1.omega_ + X2*pr2.omega_,
            Zc
        );
    }
}


template<class Specie>
inline Foam::PengRobinson<Specie> Foam::operator*
(
    const scalar s,
    const PengRobinson<Specie>& pr
)
{
    return PengRobinson<Specie>
    (
        s*static_cast<const Specie&>(pr),
        pr.Tc_,
        pr.Vc_,
        pr.Pc_,
        pr.omega_,
        pr.Zc_
    );
}


template<class Specie>
inline Foam::PengRobinson<Specie> Foam::operator==
(
    const PengRobinson<Specie>& pr1,
    const PengRobinson<Specie>& pr2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pr1)
     == static_cast<const Specie&>(pr2)
    );

    const scalar X1 = sp.W()*pr1.Y()/(pr1.W()*sp.Y());
    const scalar X2 = sp.W()*pr2.Y()/(pr2.W()*sp.Y());

    const scalar Tc = X2*pr2.Tc_ - X1*pr1.Tc_;
    const scalar Vc = X2*pr2.Vc_ - X1*pr1.Vc_;
    const scalar Zc = X2*pr2.Zc_ - X1*pr1.Zc_;

    return PengRobinson<Specie>
    (
        sp,
        Tc,
        Vc,
        RR*Zc*Tc/Vc,
        X2*pr2.omega_ - X1*pr1.omega_,
        Zc
    );
}


// ************************************************************************* //
