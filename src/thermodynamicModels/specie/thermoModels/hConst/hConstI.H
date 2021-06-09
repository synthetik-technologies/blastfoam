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

template<class EquationOfState>
inline Foam::hConst<EquationOfState>::hConst
(
    const EquationOfState& st,
    const scalar cp,
    const scalar hf,
    const scalar ft
)
:
    EquationOfState(st),
    Cp_(cp),
    Hf_(hf),
    flameT_(ft)
{
    EquationOfState::set(*this);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hConst<EquationOfState>::hConst
(
    const word& name,
    const hConst& ct
)
:
    EquationOfState(name, ct),
    Cp_(ct.Cp_),
    Hf_(ct.Hf_),
    flameT_(ct.flameT_)
{
    EquationOfState::set(*this);
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hConst<EquationOfState>>
Foam::hConst<EquationOfState>::clone() const
{
    return autoPtr<hConst<EquationOfState>>
    (
        new hConst<EquationOfState>(*this)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hConst<EquationOfState>>
Foam::hConst<EquationOfState>::New(const dictionary& dict)
{
    return autoPtr<hConst<EquationOfState>>
    (
        new hConst<EquationOfState>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::hConst<EquationOfState>::limit
(
    const scalar T
) const
{
    return max(T, 0.0);
}


template<class EquationOfState>
inline Foam::scalar Foam::hConst<EquationOfState>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Cp_ + EquationOfState::Cp(rho, e, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hConst<EquationOfState>::Hs
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Cp_*T + EquationOfState::H(rho, e, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hConst<EquationOfState>::Hf() const
{
    return Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hConst<EquationOfState>::flameT() const
{
    return flameT_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hConst<EquationOfState>::Ha
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Hs(rho, e, T) + Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hConst<EquationOfState>::S
(
    const scalar p,
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Cp_*log(T/Tstd) + EquationOfState::S(p, rho, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hConst<EquationOfState>::dGdT
(
    const scalar p,
    const scalar T
) const
{
    return 0;
}


template<class EquationOfState>
inline Foam::scalar Foam::hConst<EquationOfState>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::hConst<EquationOfState>::operator+=
(
    const hConst<EquationOfState>& ct
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(ct);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        scalar Y2 = ct.Y()/this->Y();

        Cp_ = Y1*Cp_ + Y2*ct.Cp_;
        Hf_ = Y1*Hf_ + Y2*ct.Hf_;
        flameT_ = Y1*flameT_ + Y2*ct.flameT_;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hConst<EquationOfState> Foam::operator+
(
    const hConst<EquationOfState>& ct1,
    const hConst<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    if (mag(eofs.Y()) < small)
    {
        hConst<EquationOfState> hThermo
        (
            eofs,
            ct1.Cp_,
            ct1.Hf_,
            ct1.flameT_
        );
        hThermo.EquationOfState::set(hThermo);
        return hThermo;
    }
    else
    {
        hConst<EquationOfState> hThermo
        (
            eofs,
            ct1.Y()/eofs.Y()*ct1.Cp_
          + ct2.Y()/eofs.Y()*ct2.Cp_,
            ct1.Y()/eofs.Y()*ct1.Hf_
          + ct2.Y()/eofs.Y()*ct2.Hf_,
            ct1.Y()/eofs.Y()*ct1.flameT_
          + ct2.Y()/eofs.Y()*ct2.flameT_
        );
        hThermo.EquationOfState::set(hThermo);
        return hThermo;
    }
}


template<class EquationOfState>
inline Foam::hConst<EquationOfState> Foam::operator*
(
    const scalar s,
    const hConst<EquationOfState>& ct
)
{
    hConst<EquationOfState> hThermo
    (
        s*static_cast<const EquationOfState&>(ct),
        ct.Cp_,
        ct.Hf_,
        ct.flameT_
    );
    hThermo.EquationOfState::set(hThermo);
    return hThermo;
}


template<class EquationOfState>
inline Foam::hConst<EquationOfState> Foam::operator==
(
    const hConst<EquationOfState>& ct1,
    const hConst<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
     == static_cast<const EquationOfState&>(ct2)
    );

    hConst<EquationOfState> hThermo
    (
        eofs,
        ct2.Y()/eofs.Y()*ct2.Cp_
      - ct1.Y()/eofs.Y()*ct1.Cp_,
        ct2.Y()/eofs.Y()*ct2.Hf_
      - ct1.Y()/eofs.Y()*ct1.Hf_,
        ct2.Y()/eofs.Y()*ct2.flameT_
      - ct1.Y()/eofs.Y()*ct1.flameT_
    );
    hThermo.EquationOfState::set(hThermo);
    return hThermo;
}


// ************************************************************************* //
