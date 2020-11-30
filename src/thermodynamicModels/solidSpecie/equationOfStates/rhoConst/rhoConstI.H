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
inline Foam::rhoConst<Specie>::rhoConst
(
    const Specie& sp,
    const scalar rho,
    const scalar pRef
)
:
    Specie(sp),
    rho0_(rho),
    pRef_(pRef)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::rhoConst<Specie>::rhoConst
(
    const word& name,
    const rhoConst<Specie>& ico
)
:
    Specie(name, ico),
    rho0_(ico.rho0_),
    pRef_(ico.pRef_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::rhoConst<Specie>>
Foam::rhoConst<Specie>::clone() const
{
    return autoPtr<rhoConst<Specie>>(new rhoConst<Specie>(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
Foam::scalar Foam::rhoConst<Specie>::Gamma
(
    const scalar& rho,
    const scalar& T,
    const scalar& cv
) const
{
    return 1.0;
}


template<class Specie>
Foam::scalar Foam::rhoConst<Specie>::p
(
    const scalar& rho,
    const scalar& T
) const
{
    return pRef_;
}


template<class Specie>
Foam::scalar Foam::rhoConst<Specie>::speedOfSound
(
    const scalar& p,
    const scalar& rho,
    const scalar& T,
    const scalar& cv
) const
{
    return great;
}


template<class Specie>
Foam::scalar Foam::rhoConst<Specie>::dpdv
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::rhoConst<Specie>::dpdT
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::rhoConst<Specie>::E
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::rhoConst<Specie>::Cv
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::rhoConst<Specie>::Cp
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::rhoConst<Specie>::H
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::rhoConst<Specie>::CpMCv
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::rhoConst<Specie>::S
(
    const scalar& p,
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::rhoConst<Specie>::operator+=(const rhoConst<Specie>& ico)
{
    scalar Y1 = this->Y();
    Specie::operator+=(ico);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        const scalar Y2 = ico.Y()/this->Y();

        rho0_ = Y1*rho0_ + Y2*ico.rho0_;
    }
}


template<class Specie>
inline void Foam::rhoConst<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::rhoConst<Specie> Foam::operator+
(
    const rhoConst<Specie>& ico1,
    const rhoConst<Specie>& ico2
)
{
    Specie sp
    (
        static_cast<const Specie&>(ico1)
      + static_cast<const Specie&>(ico2)
    );

    if (mag(sp.Y()) < small)
    {
        return rhoConst<Specie>
        (
            sp,
            ico1.rho0_,
            ico1.pRef_
        );
    }
    else
    {
        const scalar Y1 = ico1.Y()/sp.Y();
        const scalar Y2 = ico2.Y()/sp.Y();

        return rhoConst<Specie>
        (
            sp,
            Y1*ico1.rho0_ + Y2*ico2.rho0_,
            Y1*ico1.pRef_ + Y2*ico2.pRef_
        );
    }
}


template<class Specie>
inline Foam::rhoConst<Specie> Foam::operator*
(
    const scalar s,
    const rhoConst<Specie>& ico
)
{
    return rhoConst<Specie>
    (
        s*static_cast<const Specie&>(ico),
        ico.rho_,
        ico.pRef_
    );
}


template<class Specie>
inline Foam::rhoConst<Specie> Foam::operator==
(
    const rhoConst<Specie>& ico1,
    const rhoConst<Specie>& ico2
)
{
    Specie sp
    (
        static_cast<const Specie&>(ico1)
     == static_cast<const Specie&>(ico2)
    );

    const scalar Y1 = ico1.Y()/sp.Y();
    const scalar Y2 = ico2.Y()/sp.Y();

    return rhoConst<Specie>
    (
        sp,
        Y2*ico2.rho0_ - Y1*ico1.rho0_,
        Y2*ico2.pRef_ - Y1*ico1.pRef_
    );
}


// ************************************************************************* //
