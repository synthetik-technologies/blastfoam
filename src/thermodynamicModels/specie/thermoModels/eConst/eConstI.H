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
inline Foam::eConst<EquationOfState>::eConst
(
    const EquationOfState& st,
    const scalar cv,
    const scalar hf,
    const scalar ft
)
:
    EquationOfState(st),
    Cv_(cv),
    Hf_(hf),
    flameT_(ft)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::eConst<EquationOfState>::eConst
(
    const word& name,
    const eConst<EquationOfState>& ct
)
:
    EquationOfState(name, ct),
    Cv_(ct.Cv_),
    Hf_(ct.Hf_),
    flameT_(ct.flameT_)
{}


template<class EquationOfState>
inline Foam::autoPtr<Foam::eConst<EquationOfState>>
Foam::eConst<EquationOfState>::clone() const
{
    return autoPtr<eConst<EquationOfState>>
    (
        new eConst<EquationOfState>(*this)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::eConst<EquationOfState>>
Foam::eConst<EquationOfState>::New(const dictionary& dict)
{
    return autoPtr<eConst<EquationOfState>>
    (
        new eConst<EquationOfState>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::eConst<EquationOfState>::limit
(
    const scalar T
) const
{
    return max(T, 0.0);
}


template<class EquationOfState>
inline Foam::scalar Foam::eConst<EquationOfState>::Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Cv_ + EquationOfState::Cv(rho, e, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eConst<EquationOfState>::Es
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Cv_*T + EquationOfState::E(rho, e, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eConst<EquationOfState>::Ea
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Es(rho, e, T) + Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::eConst<EquationOfState>::Hf() const
{
    return Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::eConst<EquationOfState>::flameT() const
{
    return flameT_;
}


template<class EquationOfState>
inline Foam::scalar Foam::eConst<EquationOfState>::S
(
    const scalar p,
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Cp(rho, e, T)*log(T/Tstd) + EquationOfState::S(p, rho, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eConst<EquationOfState>::dGdT
(
    const scalar p,
    const scalar T
) const
{
    return 0;
}


template<class EquationOfState>
inline Foam::scalar Foam::eConst<EquationOfState>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    return 0;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::eConst<EquationOfState>::operator+=
(
    const eConst<EquationOfState>& ct
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(ct);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        const scalar Y2 = ct.Y()/this->Y();

        Cv_ = Y1*Cv_ + Y2*ct.Cv_;
        Hf_ = Y1*Hf_ + Y2*ct.Hf_;
        flameT_ = Y1*flameT_ + Y2*ct.flameT_;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::eConst<EquationOfState> Foam::operator+
(
    const eConst<EquationOfState>& ct1,
    const eConst<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    if (mag(eofs.Y()) < small)
    {
        return eConst<EquationOfState>
        (
            eofs,
            ct1.Cv_,
            ct1.Hf_,
            ct1.flameT_
        );
    }
    else
    {
        return eConst<EquationOfState>
        (
            eofs,
            ct1.Y()/eofs.Y()*ct1.Cv_
          + ct2.Y()/eofs.Y()*ct2.Cv_,
            ct1.Y()/eofs.Y()*ct1.Hf_
          + ct2.Y()/eofs.Y()*ct2.Hf_,
            ct1.Y()/eofs.Y()*ct1.flameT_
          + ct2.Y()/eofs.Y()*ct2.flameT_
        );
    }
}


template<class EquationOfState>
inline Foam::eConst<EquationOfState> Foam::operator*
(
    const scalar s,
    const eConst<EquationOfState>& ct
)
{
    return eConst<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(ct),
        ct.Cv_,
        ct.Hf_,
        ct.flameT_
    );
}


template<class EquationOfState>
inline Foam::eConst<EquationOfState> Foam::operator==
(
    const eConst<EquationOfState>& ct1,
    const eConst<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
     == static_cast<const EquationOfState&>(ct2)
    );

    return eConst<EquationOfState>
    (
        eofs,
        ct2.Y()/eofs.Y()*ct2.Cv_
      - ct1.Y()/eofs.Y()*ct1.Cv_,
        ct2.Y()/eofs.Y()*ct2.Hf_
      - ct1.Y()/eofs.Y()*ct1.Hf_,
        ct2.Y()/eofs.Y()*ct2.flameT_
      - ct1.Y()/eofs.Y()*ct1.flameT_
    );
}


// ************************************************************************* //
