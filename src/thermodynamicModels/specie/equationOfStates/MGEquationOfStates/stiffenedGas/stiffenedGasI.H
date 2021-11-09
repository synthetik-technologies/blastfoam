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

#include "stiffenedGas.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::stiffenedGas<Specie>::stiffenedGas
(
    const Specie& sp,
    const scalar a,
    const scalar gamma
)
:
    Specie(sp),
    a_(a),
    gamma_(gamma)
{}


template<class Specie>
Foam::scalar Foam::stiffenedGas<Specie>::Pi
(
    const scalar rho,
    const scalar e
) const
{
    return gamma_*a_;
}


template<class Specie>
Foam::scalar Foam::stiffenedGas<Specie>::delta
(
    const scalar p,
    const scalar rho,
    const scalar e
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::stiffenedGas<Specie>::GammaMG
(
    const scalar rho,
    const scalar e
) const
{
    return gamma_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::stiffenedGas<Specie>::stiffenedGas
(
    const word& name,
    const stiffenedGas<Specie>& pf
)
:
    Specie(name, pf),
    a_(pf.a_),
    gamma_(pf.gamma_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::stiffenedGas<Specie>>
Foam::stiffenedGas<Specie>::clone() const
{
    return autoPtr<stiffenedGas<Specie>>
    (
        new stiffenedGas<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::stiffenedGas<Specie>>
Foam::stiffenedGas<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<stiffenedGas<Specie>>
    (
        new stiffenedGas<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
Foam::scalar Foam::stiffenedGas<Specie>::p
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const bool limit
) const
{
    return
        limit
      ? max((gamma_ - 1.0)*rho*e - Pi(rho, e), 0.0)
      : (gamma_ - 1.0)*rho*e - Pi(rho, e);
}


template<class Specie>
Foam::scalar Foam::stiffenedGas<Specie>::Gamma
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar Cv
) const
{
    return gamma_;
}


template<class Specie>
Foam::scalar Foam::stiffenedGas<Specie>::cSqr
(
    const scalar p,
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar Cv
) const
{
    scalar h((gamma_*p + Pi(rho, e))/((gamma_ - 1.0)*max(rho, 1e-10)));
    return (h - delta(p, rho, e))*(gamma_ - 1.0);
}


template<class Specie>
Foam::scalar Foam::stiffenedGas<Specie>::dpdv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return -sqr(rho)*(gamma_ - 1.0)*e;
}


template<class Specie>
Foam::scalar Foam::stiffenedGas<Specie>::dpde
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return (gamma_ - 1.0)*rho;
}


template<class Specie>
Foam::scalar Foam::stiffenedGas<Specie>::dpdT
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
Foam::scalar Foam::stiffenedGas<Specie>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::stiffenedGas<Specie>::Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::stiffenedGas<Specie>::CpMCv
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar CpCv
) const
{
    return this->eBased_ ? (gamma_ - 1.0)*CpCv : CpCv*(gamma_ - 1.0)/gamma_;
}


template<class Specie>
Foam::scalar Foam::stiffenedGas<Specie>::E
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return a_/max(rho, minRho_);
}


template<class Specie>
Foam::scalar Foam::stiffenedGas<Specie>::H
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0;
}


template<class Specie>
Foam::scalar Foam::stiffenedGas<Specie>::S
(
    const scalar p,
    const scalar rho,
    const scalar T
) const
{
    scalar rho0((Pstd - gamma_*a_)/this->R()*Tstd);
    return this->R()*(log(rho0) - log(rho));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::stiffenedGas<Specie>::operator+=
(
    const stiffenedGas<Specie>& pf
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(pf);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = pf.Y()/this->Y();

        //- Save the original gamma
        scalar gamma1 = gamma_;

        gamma_ = 1.0/(Y1/(gamma_ - 1.0) + Y2/(pf.gamma_ - 1.0)) + 1.0;
        a_ =
            (gamma_ - 1.0)/gamma_
           *(
                Y1*a_*gamma1/(gamma1 - 1.0)
              + Y2*pf.a_*pf.gamma_/(pf.gamma_ - 1.0)
            );
    }
}


template<class Specie>
inline void Foam::stiffenedGas<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::stiffenedGas<Specie> Foam::operator+
(
    const stiffenedGas<Specie>& pf1,
    const stiffenedGas<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
      + static_cast<const Specie&>(pf2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return stiffenedGas<Specie>
        (
            sp,
            pf1.a_,
            pf1.gamma_
        );
    }
    else
    {
        const scalar Y1 = pf1.Y()/sp.Y();
        const scalar Y2 = pf2.Y()/sp.Y();

        scalar gamma =
            1.0/(Y1/(pf1.gamma_ - 1.0) + Y2/(pf2.gamma_ - 1.0)) + 1.0;

        return stiffenedGas<Specie>
        (
            sp,
            (gamma - 1.0)/gamma
           *(
                Y1*pf1.a_*pf1.gamma_/(pf1.gamma_ - 1.0)
              + Y2*pf2.a_*pf2.gamma_/(pf2.gamma_ - 1.0)
            ),
            gamma
        );
    }
}


template<class Specie>
inline Foam::stiffenedGas<Specie> Foam::operator*
(
    const scalar s,
    const stiffenedGas<Specie>& pf
)
{
    return stiffenedGas<Specie>
    (
        s*static_cast<const Specie&>(pf),
        pf.a_,
        pf.gamma_
    );
}


template<class Specie>
inline Foam::stiffenedGas<Specie> Foam::operator==
(
    const stiffenedGas<Specie>& pf1,
    const stiffenedGas<Specie>& pf2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pf1)
     == static_cast<const Specie&>(pf2)
    );

    const scalar Y1 = pf1.Y()/sp.Y();
    const scalar Y2 = pf2.Y()/sp.Y();

    // NEED TO CONFIRM
    scalar gamma =
            1.0/(Y2/(pf2.gamma_ - 1.0) - Y1/(pf1.gamma_ - 1.0)) + 1.0;
    return stiffenedGas<Specie>
    (
        sp,
        (gamma - 1.0)/gamma
       *(
            Y2*pf2.a_*pf2.gamma_/(pf2.gamma_ - 1.0)
          - Y1*pf1.a_*pf1.gamma_/(pf1.gamma_ - 1.0)
        ),
        gamma
    );
}


// ************************************************************************* //
