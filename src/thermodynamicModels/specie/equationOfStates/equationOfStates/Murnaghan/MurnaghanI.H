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
inline Foam::Murnaghan<Specie>::Murnaghan
(
    const Specie& sp,
    const scalar rho0,
    const scalar pRef,
    const scalar K0,
    const scalar kappa,
    const scalar K0Prime,
    const scalar n,
    const scalar Gamma
)
:
    Specie(sp),
    rho0_(rho0),
    pRef_(pRef),
    K0_(K0),
    kappa_(kappa),
    K0Prime_(K0Prime),
    n_(n),
    Gamma_(Gamma)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::Murnaghan<Specie>::Murnaghan
(
    const word& name,
    const Murnaghan<Specie>& pf
)
:
    Specie(name, pf),
    rho0_(pf.rho0_),
    pRef_(pf.pRef_),
    K0_(pf.K0_),
    kappa_(pf.kappa_),
    K0Prime_(pf.K0Prime_),
    n_(pf.n_),
    Gamma_(pf.Gamma_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::Murnaghan<Specie>>
Foam::Murnaghan<Specie>::clone() const
{
    return autoPtr<Murnaghan<Specie>>
    (
        new Murnaghan<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::Murnaghan<Specie>>
Foam::Murnaghan<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<Murnaghan<Specie>>
    (
        new Murnaghan<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
Foam::scalar Foam::Murnaghan<Specie>::p
(
    const scalar& rho,
    const scalar& T
) const
{
    return pRef_ + (pow(rho/rho0_, n_) - 1.0)/(kappa_*n_);
}


template<class Specie>
Foam::scalar Foam::Murnaghan<Specie>::Gamma
(
    const scalar& rho,
    const scalar& T,
    const scalar& cv
) const
{
    return Gamma_ + 1.0;
}


template<class Specie>
Foam::scalar Foam::Murnaghan<Specie>::speedOfSound
(
    const scalar& p,
    const scalar& rho,
    const scalar& T,
    const scalar& cv
) const
{
    return sqrt(pow(rho/rho0_, n_ - 1.0)/(rho0_*n_*kappa_));
}


template<class Specie>
Foam::scalar Foam::Murnaghan<Specie>::dpdv
(
    const scalar& rho,
    const scalar& T
) const
{
    return -rho*pow(rho/rho0_, n_)/kappa_;
}


template<class Specie>
Foam::scalar Foam::Murnaghan<Specie>::dpdT
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::Murnaghan<Specie>::E
(
    const scalar& rho,
    const scalar& T
) const
{
    scalar rhos(max(rho, 1e-10));
    return
       -pRef_/rhos
      + (pow(rho/rho0_, n_) - (n_ - 1.0))/(kappa_*n_*rhos*max(n_ - 1.0, small));
}


template<class Specie>
Foam::scalar Foam::Murnaghan<Specie>::Cv
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::Murnaghan<Specie>::H
(
    const scalar& rho,
    const scalar& T
) const
{
    return pow(rho/rho0_, n_ - 1.0)/(rho0_*kappa_*(n_ - 1.0));
}


template<class Specie>
Foam::scalar Foam::Murnaghan<Specie>::Cp
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::Murnaghan<Specie>::CpMCv
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::Murnaghan<Specie>::S
(
    const scalar& p,
    const scalar& rho,
    const scalar& T
) const
{
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::Murnaghan<Specie>::operator+=
(
    const Murnaghan<Specie>& pf
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(pf);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = pf.Y()/this->Y();

        rho0_ = Y1*rho0_ + Y2*pf.rho0_;
        pRef_ = Y1*pRef_ + Y2*pf.pRef_;
        K0_ = Y1*K0_ + Y2*pf.K0_;
        kappa_ = Y1*kappa_ + Y2*pf.kappa_;
        K0Prime_ = Y1*K0Prime_ + Y2*pf.K0Prime_;
        n_ = Y1*n_ + Y2*pf.n_;
        Gamma_ = Y1*Gamma_ + Y2*pf.Gamma_;
    }
}


template<class Specie>
inline void Foam::Murnaghan<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::Murnaghan<Specie> Foam::operator+
(
    const Murnaghan<Specie>& pf1,
    const Murnaghan<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
      + static_cast<const Specie&>(pf2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return Murnaghan<Specie>
        (
            sp,
            pf1.rho0_,
            pf1.pRef_,
            pf1.K0_,
            pf1.kappa_,
            pf1.K0Prime_,
            pf1.n_,
            pf1.Gamma_
        );
    }
    else
    {
        const scalar Y1 = pf1.Y()/sp.Y();
        const scalar Y2 = pf2.Y()/sp.Y();

        return Murnaghan<Specie>
        (
            sp,
            Y1*pf1.rho0_ + Y2*pf2.rho0_,
            Y1*pf1.pRef_ + Y2*pf2.pRef_,
            Y1*pf1.K0_ + Y2*pf2.K0_,
            Y1*pf1.kappa_ + Y2*pf2.kappa_,
            Y1*pf1.K0Prime_ + Y2*pf2.K0Prime_,
            Y1*pf1.n_ + Y2*pf2.n_,
            Y1*pf1.Gamma_ + Y2*pf2.Gamma_
        );
    }
}


template<class Specie>
inline Foam::Murnaghan<Specie> Foam::operator*
(
    const scalar s,
    const Murnaghan<Specie>& pf
)
{
    return Murnaghan<Specie>
    (
        s*static_cast<const Specie&>(pf),
        pf.rho0_,
        pf.pRef_,
        pf.K0_,
        pf.kappa_,
        pf.K0Prime_,
        pf.n_,
        pf.Gamma_
    );
}


template<class Specie>
inline Foam::Murnaghan<Specie> Foam::operator==
(
    const Murnaghan<Specie>& pf1,
    const Murnaghan<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
     == static_cast<const Specie&>(pf2)
    );

    const scalar Y1 = pf1.Y()/sp.Y();
    const scalar Y2 = pf2.Y()/sp.Y();

    return Murnaghan<Specie>
    (
        sp,
        Y2*pf2.rho0_ - Y1*pf1.rho0_,
        Y2*pf2.pRef_ - Y1*pf1.pRef_,
        Y2*pf2.K0_ - Y1*pf1.K0_,
        Y2*pf2.kappa_ - Y1*pf1.kappa_,
        Y2*pf2.K0Prime_ - Y1*pf1.K0Prime_,
        Y2*pf2.n_ - Y1*pf1.n_,
        Y2*pf2.Gamma_ - Y1*pf1.Gamma_
    );
}


// ************************************************************************* //
