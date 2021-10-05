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
inline Foam::eTabulatedThermo<EquationOfState>::eTabulatedThermo
(
    const EquationOfState& st,
    const scalarLookupTable1D& eTable
)
:
    EquationOfState(st),
    eTable_(eTable),
    Tlow_(min(eTable_.x())),
    Thigh_(max(eTable_.x()))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::limit
(
    const scalar T,
    const scalar Told
) const
{
    if (T/Told > 5.0)
    {
        return Told*5.0;
    }

    if
    (
        eTabulatedThermo<EquationOfState>::debug
     && (T < Tlow_ || T > Thigh_)
    )
    {
        WarningInFunction
            << "attempt to use eTabulatedThermo<EquationOfState>"
               " out of temperature range "
            << Tlow_ << " -> " << Thigh_ << ";  T = " << T
            << endl;
    }
    return min(max(T, Tlow_), Thigh_);
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return eTable_.dFdX(T) + EquationOfState::Cv(rho, e, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::Es
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return eTable_.lookup(T) + EquationOfState::E(rho, e, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::Ets
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return eTable_.lookup(T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::Ea
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Es(rho, e, T) + Hf();
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::Eta
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Ets(rho, e, T) + Hf();
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::Hf() const
{
    NotImplemented;
    return 0.0;
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::flameT() const
{
    NotImplemented;
    return 0.0;
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::S
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
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::dGdT
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::dCpdT
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    NotImplemented;
    return 0.0;
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::dCvdT
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return eTable_.d2FdX2(T);
}

// ************************************************************************* //
