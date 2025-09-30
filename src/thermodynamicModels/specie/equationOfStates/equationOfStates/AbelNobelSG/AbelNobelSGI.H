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
inline Foam::AbelNobelSG<Specie>::AbelNobelSG
(
    const Specie& sp,
    const scalar b,
    const scalar gamma,
    const scalar pInf,
    const scalar Cv
)
:
    Specie(sp),
    b_(b),
    gamma_(gamma),
    pInf_(pInf),
    Cv_(Cv)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::AbelNobelSG<Specie>::AbelNobelSG
(
    const word& name,
    const AbelNobelSG<Specie>& pf
)
:
    Specie(name, pf),
    b_(pf.b_),
    gamma_(pf.gamma_),
    pInf_(pf.pInf_),
    Cv_(pf.Cv_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::AbelNobelSG<Specie>>
Foam::AbelNobelSG<Specie>::clone() const
{
    return autoPtr<AbelNobelSG<Specie>>
    (
        new AbelNobelSG<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::AbelNobelSG<Specie>>
Foam::AbelNobelSG<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<AbelNobelSG<Specie>>
    (
        new AbelNobelSG<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
Foam::scalar Foam::AbelNobelSG<Specie>::p
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const bool limit
) const
{
    return limit ? max(pRhoT(rho, T), pCav_) : pRhoT(rho, T);
}


template<class Specie>
Foam::scalar Foam::AbelNobelSG<Specie>::pRhoT
(
    const scalar rho,
    const scalar T
) const
{
    return
        rho*(gamma_ - 1.0)*Cv_*T/(1.0 - b_*rho) - pInf_ ;
}


template<class Specie>
Foam::scalar Foam::AbelNobelSG<Specie>::rho
(
    const scalar p,
    const scalar T
) const
{
    return
        (p + pInf_)
       /((gamma_ - 1.0)*Cv_*T + b_*(p + pInf_));
}


template<class Specie>
Foam::scalar Foam::AbelNobelSG<Specie>::dpdv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return -(gamma_ - 1.0)*Cv_*T*sqr(rho/(1.0 - b_*rho));
}


template<class Specie>
Foam::scalar Foam::AbelNobelSG<Specie>::dpde
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
Foam::scalar Foam::AbelNobelSG<Specie>::dpdT
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return rho*(gamma_ - 1.0)*Cv_/(1.0 - b_*rho);
}


template<class Specie>
Foam::scalar Foam::AbelNobelSG<Specie>::Gamma
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
Foam::scalar Foam::AbelNobelSG<Specie>::cSqr
(
    const scalar p,
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar cv
) const
{
    return
        gamma_*(p + pInf_)/max(rho, 1e-10)/(1.0 - b_*rho);
//         (sqr(dpdT(rho, e, T))*T/cv - dpdv(rho, e, T))
//        /sqr(max(rho, 1e-10));
}


template<class Specie>
Foam::scalar Foam::AbelNobelSG<Specie>::E
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::AbelNobelSG<Specie>::Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::AbelNobelSG<Specie>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::AbelNobelSG<Specie>::H
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::AbelNobelSG<Specie>::CpMCv
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar CpCv
) const
{
    return this->R();
}


template<class Specie>
Foam::scalar Foam::AbelNobelSG<Specie>::Sp
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::AbelNobelSG<Specie>::Sv
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
inline void Foam::AbelNobelSG<Specie>::operator+=
(
    const AbelNobelSG<Specie>& an
)
{
    scalar X1 = this->Y()/this->W();
    Specie::operator+=(an);

    if (mag(this->Y()) > SMALL)
    {
        X1 *= this->W()/this->Y();
        const scalar X2 = this->W()*an.Y()/(an.W()*this->Y());

        b_ = X1*b_ + X2*an.b_;
    }
}


template<class Specie>
inline void Foam::AbelNobelSG<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::AbelNobelSG<Specie> Foam::operator+
(
    const AbelNobelSG<Specie>& an1,
    const AbelNobelSG<Specie>& an2
)
{
    NotImplemented;
    return an1;
//     Specie sp
//     (
//         static_cast<const Specie&>(an1)
//       + static_cast<const Specie&>(an2)
//     );
//     if (mag(sp.Y()) < SMALL)
//     {
//         return AbelNobelSG<Specie>
//         (
//             sp,
//             an1.b_
//         );
//     }
//     else
//     {
//         const scalar X1 = sp.W()*an1.Y()/(an1.W()*sp.Y());
//         const scalar X2 = sp.W()*an2.Y()/(an2.W()*sp.Y());
//
//         return AbelNobelSG<Specie>
//         (
//             sp,
//             X1*an1.b_ + X2*an2.b_
//         );
//     }
}


template<class Specie>
inline Foam::AbelNobelSG<Specie> Foam::operator*
(
    const scalar s,
    const AbelNobelSG<Specie>& pf
)
{
    NotImplemented;
    return pf;
//     return AbelNobelSG<Specie>
//     (
//         s*static_cast<const Specie&>(pf),
//         pf.b_,
//         pf.
//     );
}


template<class Specie>
inline Foam::AbelNobelSG<Specie> Foam::operator==
(
    const AbelNobelSG<Specie>& pf1,
    const AbelNobelSG<Specie>& pf2
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
//     return AbelNobelSG<Specie>
//     (
//         sp,
//         Y2*pf2.b_ - Y1*pf1.b_
//     );
}


// ************************************************************************* //
