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
inline Foam::MieGruneisen<Specie>::MieGruneisen
(
    const Specie& sp,
    const scalar rho0,
    const scalar pInf,
    const scalar c0,
    const scalar Gamma
)
:
    Specie(sp),
    rho0_(rho0),
    pInf_(pInf),
    c0_(c0),
    Gamma_(Gamma)
{}


template<class Specie>
Foam::scalar Foam::MieGruneisen<Specie>::Pi(const scalar rho) const
{
    return pInf_ + Gamma_*(rho*c0_ + pInf_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::MieGruneisen<Specie>::MieGruneisen
(
    const word& name,
    const MieGruneisen<Specie>& pf
)
:
    Specie(name, pf),
    rho0_(pf.rho0_),
    pInf_(pf.pInf_),
    c0_(pf.c0_),
    Gamma_(pf.Gamma_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::MieGruneisen<Specie>>
Foam::MieGruneisen<Specie>::clone() const
{
    return autoPtr<MieGruneisen<Specie>>
    (
        new MieGruneisen<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::MieGruneisen<Specie>>
Foam::MieGruneisen<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<MieGruneisen<Specie>>
    (
        new MieGruneisen<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
Foam::scalar Foam::MieGruneisen<Specie>::p
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const bool limit
) const
{
    scalar pMG = Gamma_*rho*e - Pi(rho);
    return limit ? max(pMG, 0.0) : pMG;
}


template<class Specie>
Foam::scalar Foam::MieGruneisen<Specie>::rho
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class Specie>
Foam::scalar Foam::MieGruneisen<Specie>::ePRho
(
    const scalar p,
    const scalar rho
) const
{
    return (p + Pi(rho))/(Gamma_*max(rho, 1e-10));
}


template<class Specie>
Foam::scalar Foam::MieGruneisen<Specie>::Gamma
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar Cv
) const
{
    return Gamma_;
}


template<class Specie>
Foam::scalar Foam::MieGruneisen<Specie>::cSqr
(
    const scalar p,
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar Cv
) const
{
    const scalar h = ((Gamma_ + 1.0)*p + Pi(rho))/(Gamma_*max(rho, 1e-10));
    return (h - c0_)*Gamma_;
}


template<class Specie>
Foam::scalar Foam::MieGruneisen<Specie>::dpdv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return -sqr(rho)*(Gamma_*(e - c0_));
}


template<class Specie>
Foam::scalar Foam::MieGruneisen<Specie>::dpde
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Gamma_*rho;
}


template<class Specie>
Foam::scalar Foam::MieGruneisen<Specie>::dpdT
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
Foam::scalar Foam::MieGruneisen<Specie>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::MieGruneisen<Specie>::Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::MieGruneisen<Specie>::CpMCv
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar CpCv
) const
{
    return this->eBased_ ? Gamma_*CpCv : CpCv*Gamma_/(Gamma_ + 1.0);
}


template<class Specie>
Foam::scalar Foam::MieGruneisen<Specie>::E
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return c0_ + pInf_/max(rho, 1e-10);
}


template<class Specie>
Foam::scalar Foam::MieGruneisen<Specie>::H
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
Foam::scalar Foam::MieGruneisen<Specie>::Sp
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::MieGruneisen<Specie>::Sv
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
inline void Foam::MieGruneisen<Specie>::operator+=
(
    const MieGruneisen<Specie>& pf
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
inline void Foam::MieGruneisen<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::MieGruneisen<Specie> Foam::operator+
(
    const MieGruneisen<Specie>& pf1,
    const MieGruneisen<Specie>& pf2
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
//         return MieGruneisen<Specie>
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
//         return MieGruneisen<Specie>
//         (
//             sp,
//             Y1*pf1.gamma_ + Y2*pf2.gamma_
//         );
//     }
}


template<class Specie>
inline Foam::MieGruneisen<Specie> Foam::operator*
(
    const scalar s,
    const MieGruneisen<Specie>& pf
)
{

    NotImplemented;
    return pf;
//     return MieGruneisen<Specie>
//     (
//         s*static_cast<const Specie&>(pf),
//         pf.gamma_
//     );
}


template<class Specie>
inline Foam::MieGruneisen<Specie> Foam::operator==
(
    const MieGruneisen<Specie>& pf1,
    const MieGruneisen<Specie>& pf2
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
//     return MieGruneisen<Specie>
//     (
//         sp,
//         Y2*pf2.gamma_ - Y1*pf1.gamma_
//     );
}


// ************************************************************************* //
