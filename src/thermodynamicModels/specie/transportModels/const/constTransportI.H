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
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::constTransport<Thermo>::constTransport
(
    const Thermo& t,
    const scalar mu,
    const scalar Pr
)
:
    Thermo(t),
    mu_(mu),
    Pr_(1.0/Pr)
{}


template<class Thermo>
inline Foam::constTransport<Thermo>::constTransport
(
    const word& name,
    const constTransport& ct
)
:
    Thermo(name, ct),
    mu_(ct.mu_),
    Pr_(ct.Pr_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::constTransport<Thermo>>
Foam::constTransport<Thermo>::clone() const
{
    return autoPtr<constTransport<Thermo>>
    (
        new constTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::constTransport<Thermo>>
Foam::constTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<constTransport<Thermo>>
    (
        new constTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::constTransport<Thermo>::mu
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return mu_;
}


template<class Thermo>
inline Foam::scalar Foam::constTransport<Thermo>::kappa
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return this->Cp(rho, e, T)*mu(rho, e, T)/Pr_;
}


template<class Thermo>
inline Foam::scalar Foam::constTransport<Thermo>::alphah
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return mu(rho, e, T)/Pr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::constTransport<Thermo>::operator+=
(
    const constTransport<Thermo>& st
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(st);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        scalar Y2 = st.Y()/this->Y();

        mu_ = Y1*mu_ + Y2*st.mu_;
        Pr_ = Y1*Pr_ + Y2*st.Pr_;
    }
}


template<class Thermo>
inline void Foam::constTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::constTransport<Thermo> Foam::operator+
(
    const constTransport<Thermo>& ct1,
    const constTransport<Thermo>& ct2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ct1) + static_cast<const Thermo&>(ct2)
    );

    if (mag(t.Y()) < small)
    {
        return constTransport<Thermo>
        (
            t,
            0,
            ct1.Pr_
        );
    }
    else
    {
        scalar Y1 = ct1.Y()/t.Y();
        scalar Y2 = ct2.Y()/t.Y();

        return constTransport<Thermo>
        (
            t,
            Y1*ct1.mu_ + Y2*ct2.mu_,
            Y1*ct1.Pr_ + Y2*ct2.Pr_
        );
    }
}


template<class Thermo>
inline Foam::constTransport<Thermo> Foam::operator*
(
    const scalar s,
    const constTransport<Thermo>& ct
)
{
    return constTransport<Thermo>
    (
        s*static_cast<const Thermo&>(ct),
        ct.mu_,
        ct.Pr_
    );
}


// ************************************************************************* //
