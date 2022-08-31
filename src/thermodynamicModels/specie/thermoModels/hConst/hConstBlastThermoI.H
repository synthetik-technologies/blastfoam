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
inline Foam::hConstThermo<EquationOfState>::hConstThermo
(
    const EquationOfState& st,
    const scalar cp,
    const scalar hf,
    const scalar Tref,
    const scalar Hsref,
    const scalar ft
)
:
    EquationOfState(st),
    Cp_(cp),
    Hf_(hf),Tref_(Tref),
    Hsref_(Hsref),
    flameT_(ft)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hConstThermo<EquationOfState>::hConstThermo
(
    const word& name,
    const hConstThermo<EquationOfState>& ct
)
:
    EquationOfState(name, ct),
    Cp_(ct.Cp_),
    Hf_(ct.Hf_),
    Tref_(ct.Tref_),
    Hsref_(ct.Hsref_),
    flameT_(ct.flameT_)
{}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hConstThermo<EquationOfState>>
Foam::hConstThermo<EquationOfState>::clone() const
{
    return autoPtr<hConstThermo<EquationOfState>>
    (
        new hConstThermo<EquationOfState>(*this)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hConstThermo<EquationOfState>>
Foam::hConstThermo<EquationOfState>::New(const dictionary& dict)
{
    return autoPtr<hConstThermo<EquationOfState>>
    (
        new hConstThermo<EquationOfState>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::limit
(
    const scalar T,
    const scalar
) const
{
    return max(T, 0);
}


template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Cp_ + EquationOfState::Cp(rho, e, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::Hs
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Cp_*(T - Tref_) + Hsref_ + EquationOfState::H(rho, e, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::Hts
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Cp_*(T - Tref_) + Hsref_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::Ha
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Hs(rho, e, T) + Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::Hta
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Hts(rho, e, T) + Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::Hf() const
{
    return Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::flameT() const
{
    return flameT_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::S
(
    const scalar p,
    const scalar T
) const
{
    return Cp_*log(T/Tstd) + EquationOfState::S(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::Gstd
(
    const scalar T
) const
{
    return Cp_*(T - Tref_) + Hsref_ + Hf() - Cp_*T*log(T/Tstd);
}


template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::dCpdT
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 0.0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::hConstThermo<EquationOfState>::operator+=
(
    const hConstThermo<EquationOfState>& ct
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(ct);

    if (mag(this->Y()) > small)
    {
        if
        (
            hConstThermo<EquationOfState>::debug
         && notEqual(Tref_, ct.Tref_)
        )
        {
            FatalErrorInFunction
                << "Tref " << Tref_ << " for "
                << (this->name().size() ? this->name() : "others")
                << " != " << ct.Tref_ << " for "
                << (ct.name().size() ? ct.name() : "others")
                << exit(FatalError);
        }

        Y1 /= this->Y();
        const scalar Y2 = ct.Y()/this->Y();

        Cp_ = Y1*Cp_ + Y2*ct.Cp_;
        Hf_ = Y1*Hf_ + Y2*ct.Hf_;
        Hsref_ = Y1*Hsref_ + Y2*ct.Hsref_;
        flameT_ = Y1*flameT_ + Y2*ct.flameT_;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hConstThermo<EquationOfState> Foam::operator+
(
    const hConstThermo<EquationOfState>& ct1,
    const hConstThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    if (mag(eofs.Y()) < small)
    {
        return hConstThermo<EquationOfState>
        (
            eofs,
            ct1.Cp_,
            ct1.Hf_,
            ct1.Tref_,
            ct1.Hsref_,
            ct1.flameT_
        );
    }
    else
    {
        if
        (
            hConstThermo<EquationOfState>::debug
         && notEqual(ct1.Tref_, ct2.Tref_)
        )
        {
            FatalErrorInFunction
                << "Tref " << ct1.Tref_ << " for "
                << (ct1.name().size() ? ct1.name() : "others")
                << " != " << ct2.Tref_ << " for "
                << (ct2.name().size() ? ct2.name() : "others")
                << exit(FatalError);
        }

        return hConstThermo<EquationOfState>
        (
            eofs,
            ct1.Y()/eofs.Y()*ct1.Cp_
          + ct2.Y()/eofs.Y()*ct2.Cp_,
            ct1.Y()/eofs.Y()*ct1.Hf_
          + ct2.Y()/eofs.Y()*ct2.Hf_,
            ct1.Tref_,
            ct1.Y()/eofs.Y()*ct1.Hsref_
          + ct2.Y()/eofs.Y()*ct2.Hsref_,
            ct1.Y()/eofs.Y()*ct1.flameT_
          + ct2.Y()/eofs.Y()*ct2.flameT_
        );
    }
}


template<class EquationOfState>
inline Foam::hConstThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const hConstThermo<EquationOfState>& ct
)
{
    return hConstThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(ct),
        ct.Cp_,
        ct.Hf_,
        ct.Tref_,
        ct.Hsref_,
        ct.flameT_
    );
}


template<class EquationOfState>
inline Foam::hConstThermo<EquationOfState> Foam::operator==
(
    const hConstThermo<EquationOfState>& ct1,
    const hConstThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
     == static_cast<const EquationOfState&>(ct2)
    );

    if
    (
        hConstThermo<EquationOfState>::debug
     && notEqual(ct1.Tref_, ct2.Tref_)
    )
    {
        FatalErrorInFunction
            << "Tref " << ct1.Tref_ << " for "
            << (ct1.name().size() ? ct1.name() : "others")
            << " != " << ct2.Tref_ << " for "
            << (ct2.name().size() ? ct2.name() : "others")
            << exit(FatalError);
    }

    return hConstThermo<EquationOfState>
    (
        eofs,
        ct2.Y()/eofs.Y()*ct2.Cp_
      - ct1.Y()/eofs.Y()*ct1.Cp_,
        ct2.Y()/eofs.Y()*ct2.Hf_
      - ct1.Y()/eofs.Y()*ct1.Hf_,
        ct1.Tref_,
        ct2.Y()/eofs.Y()*ct2.Hsref_
      - ct1.Y()/eofs.Y()*ct1.Hsref_,
        ct2.Y()/eofs.Y()*ct2.flameT_
      - ct1.Y()/eofs.Y()*ct1.flameT_
    );
}


// ************************************************************************* //
