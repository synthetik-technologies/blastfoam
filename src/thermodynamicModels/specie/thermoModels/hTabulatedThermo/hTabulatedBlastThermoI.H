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
inline Foam::hTabulatedThermo<EquationOfState>::hTabulatedThermo
(
    const EquationOfState& st,
    const scalarLookupTable2D& hTable,
    const scalar Hf
)
:
    EquationOfState(st),
    hTable_(hTable),
    Tlow_(min(hTable_.x())),
    Thigh_(max(hTable_.x())),
    Hf_(Hf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::hTabulatedThermo<EquationOfState>::limit
(
    const scalar T,
    const scalar Told
) const
{
    if (T/Told > 5.0)
    {
        return Told*5.0;
    }

    if (hTabulatedThermo<EquationOfState>::debug && (T < Tlow_ || T > Thigh_))
    {
        WarningInFunction
            << "attempt to use hTabulatedThermo<EquationOfState>"
               " out of temperature range "
            << Tlow_ << " -> " << Thigh_ << ";  T = " << T
            << endl;
    }
    return min(max(T, Tlow_), Thigh_);
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabulatedThermo<EquationOfState>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return hTable_.dFdY(rho, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabulatedThermo<EquationOfState>::Hs
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return hTable_.lookup(rho, T) + EquationOfState::H(rho, e, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabulatedThermo<EquationOfState>::Hts
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return hTable_.lookup(rho, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabulatedThermo<EquationOfState>::Ha
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Hs(rho, e, T) + Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabulatedThermo<EquationOfState>::Hta
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Hts(rho, e, T) + Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabulatedThermo<EquationOfState>::Hf() const
{
    return Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabulatedThermo<EquationOfState>::flameT() const
{
    NotImplemented;
    return 0.0;
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabulatedThermo<EquationOfState>::S
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabulatedThermo<EquationOfState>::Gstd
(
    const scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabulatedThermo<EquationOfState>::dCpdT
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return hTable_.d2FdX2(rho, T);
}


// ************************************************************************* //
