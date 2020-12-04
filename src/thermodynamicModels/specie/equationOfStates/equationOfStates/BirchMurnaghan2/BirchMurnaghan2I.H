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
inline Foam::BirchMurnaghan2<Specie>::BirchMurnaghan2
(
    const Specie& sp,
    const scalar rho0,
    const scalar pRef,
    const scalar K0,
    const scalar K0Prime,
    const scalar Gamma
)
:
    Specie(sp),
    rho0_(rho0),
    pRef_(pRef),
    K0_(K0),
    K0Prime_(K0Prime),
    Gamma_(Gamma)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::BirchMurnaghan2<Specie>::BirchMurnaghan2
(
    const word& name,
    const BirchMurnaghan2<Specie>& pf
)
:
    Specie(name, pf),
    rho0_(pf.rho0_),
    pRef_(pf.pRef_),
    K0_(pf.K0_),
    K0Prime_(pf.K0Prime_),
    Gamma_(pf.Gamma_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::BirchMurnaghan2<Specie>>
Foam::BirchMurnaghan2<Specie>::clone() const
{
    return autoPtr<BirchMurnaghan2<Specie>>
    (
        new BirchMurnaghan2<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::BirchMurnaghan2<Specie>>
Foam::BirchMurnaghan2<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<BirchMurnaghan2<Specie>>
    (
        new BirchMurnaghan2<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
Foam::scalar Foam::BirchMurnaghan2<Specie>::p
(
    const scalar& rho,
    const scalar& T
) const
{
    return  pRef_ + 1.5*K0_*(pow(rho/rho0_, 7.0/3.0) - pow(rho/rho0_, 5.0/3.0));
}


template<class Specie>
Foam::scalar Foam::BirchMurnaghan2<Specie>::Gamma
(
    const scalar& rho,
    const scalar& T,
    const scalar& cv
) const
{
    return Gamma_ + 1.0;
}


template<class Specie>
Foam::scalar Foam::BirchMurnaghan2<Specie>::cSqr
(
    const scalar& p,
    const scalar& rho,
    const scalar& T,
    const scalar& cv
) const
{
    return
        K0_
       *pow(rho/rho0_, 5.0/3.0)
       *(7.0*pow(rho/rho0_, 2.0/3.0) - 5.0)
       /(2.0*max(rho, 1e-10));
}


template<class Specie>
Foam::scalar Foam::BirchMurnaghan2<Specie>::dpdv
(
    const scalar& rho,
    const scalar& T
) const
{
    return
        0.5*K0_*rho
       *pow(rho/rho0_, 5.0/3.0)
       *(5.0 - 7.0*pow(rho/rho0_, 2.0/3.0));
}


template<class Specie>
Foam::scalar Foam::BirchMurnaghan2<Specie>::dpdT
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0;
}


template<class Specie>
Foam::scalar Foam::BirchMurnaghan2<Specie>::E
(
    const scalar& rho,
    const scalar& T
) const
{
    return
        (
            9.0/8.0*K0_*pow(rho/rho0_, 5.0/3.0)*(pow(rho/rho0_, 2.0/3.0) - 2.0)
          - pRef_
        )/max(rho, 1e-10);
}


template<class Specie>
Foam::scalar Foam::BirchMurnaghan2<Specie>::Cv
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0;
}


template<class Specie>
Foam::scalar Foam::BirchMurnaghan2<Specie>::H
(
    const scalar& rho,
    const scalar& T
) const
{
    return 3.0/8.0*K0_*pow(rho/rho0_, 2.0/3.0)*(7.0*pow(rho/rho0_, 2.0/3.0) - 10.0);
}


template<class Specie>
Foam::scalar Foam::BirchMurnaghan2<Specie>::Cp
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::BirchMurnaghan2<Specie>::CpMCv
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::BirchMurnaghan2<Specie>::S
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
inline void Foam::BirchMurnaghan2<Specie>::operator+=
(
    const BirchMurnaghan2<Specie>& pf
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
        K0Prime_ = Y1*K0Prime_ + Y2*pf.K0Prime_;
        Gamma_ = Y1*Gamma_ + Y2*pf.Gamma_;
    }
}


template<class Specie>
inline void Foam::BirchMurnaghan2<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::BirchMurnaghan2<Specie> Foam::operator+
(
    const BirchMurnaghan2<Specie>& pf1,
    const BirchMurnaghan2<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
      + static_cast<const Specie&>(pf2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return BirchMurnaghan2<Specie>
        (
            sp,
            pf1.rho0_,
            pf1.pRef_,
            pf1.K0_,
            pf1.K0Prime_,
            pf1.Gamma_
        );
    }
    else
    {
        const scalar Y1 = pf1.Y()/sp.Y();
        const scalar Y2 = pf2.Y()/sp.Y();

        return BirchMurnaghan2<Specie>
        (
            sp,
            Y1*pf1.rho0_ + Y2*pf2.rho0_,
            Y1*pf1.pRef_ + Y2*pf2.pRef_,
            Y1*pf1.K0_ + Y2*pf2.K0_,
            Y1*pf1.K0Prime_ + Y2*pf2.K0Prime_,
            Y1*pf1.Gamma_ + Y2*pf2.Gamma_
        );
    }
}


template<class Specie>
inline Foam::BirchMurnaghan2<Specie> Foam::operator*
(
    const scalar s,
    const BirchMurnaghan2<Specie>& pf
)
{
    return BirchMurnaghan2<Specie>
    (
        s*static_cast<const Specie&>(pf),
        pf.rho0_,
        pf.pRef_,
        pf.K0_,
        pf.K0Prime_,
        pf.Gamma_
    );
}


template<class Specie>
inline Foam::BirchMurnaghan2<Specie> Foam::operator==
(
    const BirchMurnaghan2<Specie>& pf1,
    const BirchMurnaghan2<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
     == static_cast<const Specie&>(pf2)
    );

    const scalar Y1 = pf1.Y()/sp.Y();
    const scalar Y2 = pf2.Y()/sp.Y();

    return BirchMurnaghan2<Specie>
    (
        sp,
        Y2*pf2.rho0_ - Y1*pf1.rho0_,
        Y2*pf2.pRef_ - Y1*pf1.pRef_,
        Y2*pf2.K0_ - Y1*pf1.K0_,
        Y2*pf2.K0Prime_ - Y1*pf1.K0Prime_,
        Y2*pf2.Gamma_ - Y1*pf1.Gamma_
    );
}


// ************************************************************************* //
