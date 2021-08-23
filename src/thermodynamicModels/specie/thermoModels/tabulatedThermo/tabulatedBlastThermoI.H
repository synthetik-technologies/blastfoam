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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::tabulatedThermo<EquationOfState>::tabulatedThermo
(
    const tabulatedThermo& tt
)
:
    EquationOfState(tt),
    CpCvModel(CpCvModel::InternalEnergy),
    eTable_(tt.eTable_),
    Tlow_(tt.Tlow_),
    Thigh_(tt.Thigh_)
{}


template<class EquationOfState>
inline Foam::scalar Foam::tabulatedThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    if (tabulatedThermo<EquationOfState>::debug && (T < Tlow_ || T > Thigh_))
    {
        WarningInFunction
            << "attempt to use tabulatedThermo<EquationOfState>"
               " out of temperature range "
            << Tlow_ << " -> " << Thigh_ << ";  T = " << T
            << endl;

        return min(max(T, Tlow_), Thigh_);
    }
    else
    {
        return T;
    }
}


template<class EquationOfState>
inline Foam::scalar Foam::tabulatedThermo<EquationOfState>::Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return eTable_.dFdX(T) + EquationOfState::Cv(rho, e, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::tabulatedThermo<EquationOfState>::Es
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return eTable_.lookup(T);
}


template<class EquationOfState>
inline Foam::scalar Foam::tabulatedThermo<EquationOfState>::Ets
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return eTable_.lookup(T);
}


template<class EquationOfState>
inline Foam::scalar Foam::tabulatedThermo<EquationOfState>::Ea
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Es(rho, e, T) + Hf();
}


template<class EquationOfState>
inline Foam::scalar Foam::tabulatedThermo<EquationOfState>::Eta
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Ets(rho, e, T) + Hf();
}


template<class EquationOfState>
inline Foam::scalar Foam::tabulatedThermo<EquationOfState>::Hf() const
{
    NotImplemented;
    return 0.0;
}


template<class EquationOfState>
inline Foam::scalar Foam::tabulatedThermo<EquationOfState>::flameT() const
{
    NotImplemented;
    return 0.0;
}


template<class EquationOfState>
inline Foam::scalar Foam::tabulatedThermo<EquationOfState>::S
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
inline Foam::scalar Foam::tabulatedThermo<EquationOfState>::dGdT
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class EquationOfState>
inline Foam::scalar Foam::tabulatedThermo<EquationOfState>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return 0;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::tabulatedThermo<EquationOfState>::operator=
(
    const tabulatedThermo<EquationOfState>& ct
)
{
    NotImplemented;
    return;

//     scalar Y1 = this->Y();
//
//     EquationOfState::operator+=(ct);
//
//     if (mag(this->Y()) > small)
//     {
//         Y1 /= this->Y();
//         const scalar Y2 = ct.Y()/this->Y();
//
//         Cv_ = Y1*Cv_ + Y2*ct.Cv_;
//         Hf_ = Y1*Hf_ + Y2*ct.Hf_;
//     }
}

template<class EquationOfState>
inline void Foam::tabulatedThermo<EquationOfState>::operator+=
(
    const tabulatedThermo<EquationOfState>& ct
)
{
    NotImplemented;
    return;

//     scalar Y1 = this->Y();
//
//     EquationOfState::operator+=(ct);
//
//     if (mag(this->Y()) > small)
//     {
//         Y1 /= this->Y();
//         const scalar Y2 = ct.Y()/this->Y();
//
//         Cv_ = Y1*Cv_ + Y2*ct.Cv_;
//         Hf_ = Y1*Hf_ + Y2*ct.Hf_;
//     }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::tabulatedThermo<EquationOfState> Foam::operator+
(
    const tabulatedThermo<EquationOfState>& ct1,
    const tabulatedThermo<EquationOfState>& ct2
)
{
    NotImplemented;
    return ct1;

//     EquationOfState eofs
//     (
//         static_cast<const EquationOfState&>(ct1)
//       + static_cast<const EquationOfState&>(ct2)
//     );
//
//     if (mag(eofs.Y()) < small)
//     {
//         return tabulatedThermo<EquationOfState>
//         (
//             eofs,
//             ct1.Cv_,
//             ct1.Hf_
//         );
//     }
//     else
//     {
//         return tabulatedThermo<EquationOfState>
//         (
//             eofs,
//             ct1.Y()/eofs.Y()*ct1.Cv_
//           + ct2.Y()/eofs.Y()*ct2.Cv_,
//             ct1.Y()/eofs.Y()*ct1.Hf_
//           + ct2.Y()/eofs.Y()*ct2.Hf_
//         );
//     }
}


template<class EquationOfState>
inline Foam::tabulatedThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const tabulatedThermo<EquationOfState>& ct
)
{
    NotImplemented;
    return ct;

//     return tabulatedThermo<EquationOfState>
//     (
//         s*static_cast<const EquationOfState&>(ct),
//         ct.Cv_,
//         ct.Hf_
//     );
}


template<class EquationOfState>
inline Foam::tabulatedThermo<EquationOfState> Foam::operator==
(
    const tabulatedThermo<EquationOfState>& ct1,
    const tabulatedThermo<EquationOfState>& ct2
)
{
    NotImplemented;
    return ct1;

//     EquationOfState eofs
//     (
//         static_cast<const EquationOfState&>(ct1)
//      == static_cast<const EquationOfState&>(ct2)
//     );
//
//     return tabulatedThermo<EquationOfState>
//     (
//         eofs,
//         ct2.Y()/eofs.Y()*ct2.Cv_
//       - ct1.Y()/eofs.Y()*ct1.Cv_,
//         ct2.Y()/eofs.Y()*ct2.Hf_
//       - ct1.Y()/eofs.Y()*ct1.Hf_
//     );
}


// ************************************************************************* //
