/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
inline Foam::eMeltingThermo<EquationOfState>::eMeltingThermo
(
    const EquationOfState& st,
    const scalar Cv,
    const scalar CvLiquid,
    const scalar Tm,
    const scalar Lm,
    const scalar dTm,
    const scalar Hf,
    const scalar Tref,
    const scalar Esref
)
:
    EquationOfState(st),
    Cv_(Cv),
    CvLiquid_(CvLiquid),
    Tm_(Tm),
    Lm_(Lm),
    dTm_(dTm),
    Hf_(Hf),
    Tref_(Tref),
    Esref_(Esref)
{}


template<class EquationOfState>
inline Foam::scalar Foam::eMeltingThermo<EquationOfState>::T1() const
{
    return Tm_ - 0.5*dTm_;
}


template<class EquationOfState>
inline Foam::scalar Foam::eMeltingThermo<EquationOfState>::T2() const
{
    return Tm_ + 0.5*dTm_;
}


template<class EquationOfState>
inline Foam::scalar Foam::eMeltingThermo<EquationOfState>::dEt
(
    const scalar T
) const
{
    return
        Cv_*(min(T, T1()) - Tref_)
      + (Cv_ + Lm_/dTm_)*min(max(T - T1(), 0.0), dTm_)
      + CvLiquid_*max(T - T2(), 0.0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::eMeltingThermo<EquationOfState>::eMeltingThermo
(
    const word& name,
    const eMeltingThermo<EquationOfState>& mt
)
:
    EquationOfState(name, mt),
    Cv_(mt.Cv_),
    CvLiquid_(mt.CvLiquid_),
    Tm_(mt.Tm_),
    Lm_(mt.Lm_),
    dTm_(mt.dTm_),
    Hf_(mt.Hf_),
    Tref_(mt.Tref_),
    Esref_(mt.Esref_)
{}


template<class EquationOfState>
inline Foam::autoPtr<Foam::eMeltingThermo<EquationOfState>>
Foam::eMeltingThermo<EquationOfState>::clone() const
{
    return autoPtr<eMeltingThermo<EquationOfState>>
    (
        new eMeltingThermo<EquationOfState>(*this)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::eMeltingThermo<EquationOfState>>
Foam::eMeltingThermo<EquationOfState>::New(const dictionary& dict)
{
    return autoPtr<eMeltingThermo<EquationOfState>>
    (
        new eMeltingThermo<EquationOfState>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::eMeltingThermo<EquationOfState>::limit
(
    const scalar T,
    const scalar
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::eMeltingThermo<EquationOfState>::Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Cvt(T) + EquationOfState::Cv(rho, e, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eMeltingThermo<EquationOfState>::Es
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Est(T) + EquationOfState::E(rho, e, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eMeltingThermo<EquationOfState>::Ea
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Es(rho, e, T) + Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::eMeltingThermo<EquationOfState>::TEsEstimate
(
    const scalar rho,
    const scalar es,
    const scalar T0
) const
{
    // Thermal energy relative to Tref
    const scalar de = es - EquationOfState::E(rho, es, T0) - Esref_;

    const scalar de1 = Cv_*(T1() - Tref_);
    const scalar de2 = de1 + Cv_*dTm_ + Lm_;

    if (de <= de1)
    {
        return Tref_ + de/Cv_;
    }
    else if (de <= de2)
    {
        return T1() + (de - de1)/(Cv_ + Lm_/dTm_);
    }
    return T2() + (de - de2)/CvLiquid_;
}


template<class EquationOfState>
inline Foam::scalar Foam::eMeltingThermo<EquationOfState>::Cvt
(
    const scalar T
) const
{
    if (T <= T1())
    {
        return Cv_;
    }
    else if (T <= T2())
    {
        return Cv_ + Lm_/dTm_;
    }
    return CvLiquid_;
}


template<class EquationOfState>
inline Foam::scalar Foam::eMeltingThermo<EquationOfState>::Est
(
    const scalar T
) const
{
    return Esref_ + dEt(T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eMeltingThermo<EquationOfState>::Eat
(
    const scalar T
) const
{
    return Est(T) + Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::eMeltingThermo<EquationOfState>::Hf() const
{
    return Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::eMeltingThermo<EquationOfState>::S
(
    const scalar p,
    const scalar T
) const
{
    // Integral of Cv/T from Tstd, piecewise
    const scalar Ts = max(T, small);
    const scalar Tsol = min(Ts, T1());
    const scalar Tmelt = min(max(Ts, T1()), T2());

    return
        Cv_*log(Tsol/Tstd)
      + (Cv_ + Lm_/dTm_)*log(Tmelt/T1())
      + CvLiquid_*log(max(Ts, T2())/T2())
      + EquationOfState::Sv(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eMeltingThermo<EquationOfState>::Gstd
(
    const scalar T
) const
{
    return
        Est(T) + Hf_ + Pstd/EquationOfState::rho(Pstd, T)
      - S(Pstd, T)*T;
}


template<class EquationOfState>
inline Foam::scalar Foam::eMeltingThermo<EquationOfState>::dCpdT
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    NotImplemented;
    return 0.0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::eMeltingThermo<EquationOfState>::operator+=
(
    const eMeltingThermo<EquationOfState>& mt
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(mt);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        const scalar Y2 = mt.Y()/this->Y();

        Cv_ = Y1*Cv_ + Y2*mt.Cv_;
        CvLiquid_ = Y1*CvLiquid_ + Y2*mt.CvLiquid_;
        Lm_ = Y1*Lm_ + Y2*mt.Lm_;
        Hf_ = Y1*Hf_ + Y2*mt.Hf_;
        Esref_ = Y1*Esref_ + Y2*mt.Esref_;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::eMeltingThermo<EquationOfState> Foam::operator+
(
    const eMeltingThermo<EquationOfState>& mt1,
    const eMeltingThermo<EquationOfState>& mt2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(mt1)
      + static_cast<const EquationOfState&>(mt2)
    );

    if (mag(eofs.Y()) < small)
    {
        return eMeltingThermo<EquationOfState>
        (
            eofs,
            mt1.Cv_,
            mt1.CvLiquid_,
            mt1.Tm_,
            mt1.Lm_,
            mt1.dTm_,
            mt1.Hf_,
            mt1.Tref_,
            mt1.Esref_
        );
    }

    const scalar Y1 = mt1.Y()/eofs.Y();
    const scalar Y2 = mt2.Y()/eofs.Y();

    return eMeltingThermo<EquationOfState>
    (
        eofs,
        Y1*mt1.Cv_ + Y2*mt2.Cv_,
        Y1*mt1.CvLiquid_ + Y2*mt2.CvLiquid_,
        mt1.Tm_,
        Y1*mt1.Lm_ + Y2*mt2.Lm_,
        mt1.dTm_,
        Y1*mt1.Hf_ + Y2*mt2.Hf_,
        mt1.Tref_,
        Y1*mt1.Esref_ + Y2*mt2.Esref_
    );
}


template<class EquationOfState>
inline Foam::eMeltingThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const eMeltingThermo<EquationOfState>& mt
)
{
    return eMeltingThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(mt),
        mt.Cv_,
        mt.CvLiquid_,
        mt.Tm_,
        mt.Lm_,
        mt.dTm_,
        mt.Hf_,
        mt.Tref_,
        mt.Esref_
    );
}


template<class EquationOfState>
inline Foam::eMeltingThermo<EquationOfState> Foam::operator==
(
    const eMeltingThermo<EquationOfState>& mt1,
    const eMeltingThermo<EquationOfState>& mt2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(mt1)
     == static_cast<const EquationOfState&>(mt2)
    );

    const scalar Y1 = mt1.Y()/eofs.Y();
    const scalar Y2 = mt2.Y()/eofs.Y();

    return eMeltingThermo<EquationOfState>
    (
        eofs,
        Y2*mt2.Cv_ - Y1*mt1.Cv_,
        Y2*mt2.CvLiquid_ - Y1*mt1.CvLiquid_,
        mt1.Tm_,
        Y2*mt2.Lm_ - Y1*mt1.Lm_,
        mt1.dTm_,
        Y2*mt2.Hf_ - Y1*mt1.Hf_,
        mt1.Tref_,
        Y2*mt2.Esref_ - Y1*mt1.Esref_
    );
}


// ************************************************************************* //
