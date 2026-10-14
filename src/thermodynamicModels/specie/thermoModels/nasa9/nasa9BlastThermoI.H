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
inline Foam::label Foam::nasa9Thermo<EquationOfState>::range
(
    const scalar T
) const
{
    const label n = coeffs_.size();
    for (label i = 0; i < n - 1; i++)
    {
        if (T < Tswitch_[i + 1])
        {
            return i;
        }
    }
    return n - 1;
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9Thermo<EquationOfState>::Tfit
(
    const scalar T
) const
{
    return min(max(T, Tswitch_.first()), Tswitch_.last());
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9Thermo<EquationOfState>::CpFit
(
    const scalar Tf
) const
{
    const coeffArray& a = coeffs_[range(Tf)];
    const scalar rT = 1.0/Tf;
    return
        rT*(rT*a[0] + a[1])
      + a[2] + Tf*(a[3] + Tf*(a[4] + Tf*(a[5] + Tf*a[6])));
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9Thermo<EquationOfState>::HaFit
(
    const scalar Tf
) const
{
    const coeffArray& a = coeffs_[range(Tf)];
    return
      - a[0]/Tf + a[1]*log(Tf) + a[7]
      + Tf*(a[2] + Tf*(a[3]/2.0 + Tf*(a[4]/3.0 + Tf*(a[5]/4.0 + Tf*a[6]/5.0))));
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9Thermo<EquationOfState>::SFit
(
    const scalar Tf
) const
{
    const coeffArray& a = coeffs_[range(Tf)];
    const scalar rT = 1.0/Tf;
    return
      - rT*(0.5*a[0]*rT + a[1]) + a[2]*log(Tf) + a[8]
      + Tf*(a[3] + Tf*(a[4]/2.0 + Tf*(a[5]/3.0 + Tf*a[6]/4.0)));
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9Thermo<EquationOfState>::Sstd
(
    const scalar T
) const
{
    const scalar Tf = Tfit(T);
    return SFit(Tf) + CpFit(Tf)*log(T/Tf);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::nasa9Thermo<EquationOfState>::nasa9Thermo
(
    const EquationOfState& st,
    const List<scalar>& Tswitch,
    const List<coeffArray>& coeffs,
    const bool convertCoeffs
)
:
    EquationOfState(st),
    Tswitch_(Tswitch),
    coeffs_(coeffs)
{
    if (convertCoeffs)
    {
        forAll(coeffs_, i)
        {
            for (label j = 0; j < nCoeffs_; j++)
            {
                coeffs_[i][j] *= this->R();
            }
        }
    }
}


template<class EquationOfState>
inline Foam::nasa9Thermo<EquationOfState>::nasa9Thermo
(
    const word& name,
    const nasa9Thermo& nt
)
:
    EquationOfState(name, nt),
    Tswitch_(nt.Tswitch_),
    coeffs_(nt.coeffs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::nasa9Thermo<EquationOfState>::limit
(
    const scalar T,
    const scalar Told
) const
{
    if (T < Tswitch_.first() || T > Tswitch_.last())
    {
        if (nasa9Thermo<EquationOfState>::debug)
        {
            WarningInFunction
                << "attempt to use nasa9Thermo<EquationOfState>"
                << " out of temperature range "
                << Tswitch_.first() << " -> " << Tswitch_.last()
                << ";  T = " << T << endl;
        }

        // Above Thigh the properties are extrapolated with constant Cp; the
        // upper bound keeps the energy -> temperature iteration bounded in
        // non-physical states
        return min(max(T, Tswitch_.first()), 10*Tswitch_.last());
    }
    else
    {
        return T;
    }
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9Thermo<EquationOfState>::Tlow() const
{
    return Tswitch_.first();
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9Thermo<EquationOfState>::Thigh() const
{
    return Tswitch_.last();
}


template<class EquationOfState>
inline const Foam::List<Foam::scalar>&
Foam::nasa9Thermo<EquationOfState>::Tswitch() const
{
    return Tswitch_;
}


template<class EquationOfState>
inline const Foam::List
<
    typename Foam::nasa9Thermo<EquationOfState>::coeffArray
>&
Foam::nasa9Thermo<EquationOfState>::coeffs() const
{
    return coeffs_;
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9Thermo<EquationOfState>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Cpt(T) + EquationOfState::Cp(rho, e, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9Thermo<EquationOfState>::Ha
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Hat(T) + EquationOfState::H(rho, e, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9Thermo<EquationOfState>::Hs
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Ha(rho, e, T) - Hf();
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9Thermo<EquationOfState>::Cpt
(
    const scalar T
) const
{
    return CpFit(Tfit(T));
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9Thermo<EquationOfState>::Hat
(
    const scalar T
) const
{
    const scalar Tf = Tfit(T);
    return HaFit(Tf) + CpFit(Tf)*(T - Tf);
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9Thermo<EquationOfState>::Hst
(
    const scalar T
) const
{
    return Hat(T) - Hf();
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9Thermo<EquationOfState>::Hf() const
{
    return HaFit(Tstd);
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9Thermo<EquationOfState>::S
(
    const scalar p,
    const scalar T
) const
{
    return Sstd(T) + EquationOfState::Sp(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9Thermo<EquationOfState>::Gstd
(
    const scalar T
) const
{
    return Hat(T) - T*Sstd(T);
}


template<class EquationOfState>
inline Foam::scalar Foam::nasa9Thermo<EquationOfState>::dCpdT
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (T < Tswitch_.first() || T > Tswitch_.last())
    {
        return 0;
    }
    const coeffArray& a = coeffs_[range(T)];
    const scalar rT = 1.0/T;
    return
      - rT*rT*(2*a[0]*rT + a[1])
      + a[3] + T*(2*a[4] + T*(3*a[5] + T*4*a[6]));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::nasa9Thermo<EquationOfState>::operator+=
(
    const nasa9Thermo<EquationOfState>& nt
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(nt);

    if (mag(this->Y()) > small)
    {
        checkRanges(*this, nt);

        Y1 /= this->Y();
        const scalar Y2 = nt.Y()/this->Y();

        forAll(coeffs_, i)
        {
            for (label j = 0; j < nCoeffs_; j++)
            {
                coeffs_[i][j] = Y1*coeffs_[i][j] + Y2*nt.coeffs_[i][j];
            }
        }
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::nasa9Thermo<EquationOfState> Foam::operator+
(
    const nasa9Thermo<EquationOfState>& nt1,
    const nasa9Thermo<EquationOfState>& nt2
)
{
    EquationOfState eofs = nt1;
    eofs += nt2;

    if (mag(eofs.Y()) < small)
    {
        return nasa9Thermo<EquationOfState>
        (
            eofs,
            nt1.Tswitch_,
            nt1.coeffs_
        );
    }
    else
    {
        nasa9Thermo<EquationOfState>::checkRanges(nt1, nt2);

        const scalar Y1 = nt1.Y()/eofs.Y();
        const scalar Y2 = nt2.Y()/eofs.Y();

        List<typename nasa9Thermo<EquationOfState>::coeffArray>
            coeffs(nt1.coeffs_.size());

        forAll(coeffs, i)
        {
            for (label j = 0; j < nasa9Thermo<EquationOfState>::nCoeffs_; j++)
            {
                coeffs[i][j] = Y1*nt1.coeffs_[i][j] + Y2*nt2.coeffs_[i][j];
            }
        }

        return nasa9Thermo<EquationOfState>(eofs, nt1.Tswitch_, coeffs);
    }
}


template<class EquationOfState>
inline Foam::nasa9Thermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const nasa9Thermo<EquationOfState>& nt
)
{
    return nasa9Thermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(nt),
        nt.Tswitch_,
        nt.coeffs_
    );
}


template<class EquationOfState>
inline Foam::nasa9Thermo<EquationOfState> Foam::operator==
(
    const nasa9Thermo<EquationOfState>& nt1,
    const nasa9Thermo<EquationOfState>& nt2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(nt1)
     == static_cast<const EquationOfState&>(nt2)
    );

    nasa9Thermo<EquationOfState>::checkRanges(nt1, nt2);

    const scalar Y1 = nt2.Y()/eofs.Y();
    const scalar Y2 = nt1.Y()/eofs.Y();

    List<typename nasa9Thermo<EquationOfState>::coeffArray>
        coeffs(nt1.coeffs_.size());

    forAll(coeffs, i)
    {
        for (label j = 0; j < nasa9Thermo<EquationOfState>::nCoeffs_; j++)
        {
            coeffs[i][j] = Y1*nt2.coeffs_[i][j] - Y2*nt1.coeffs_[i][j];
        }
    }

    return nasa9Thermo<EquationOfState>(eofs, nt2.Tswitch_, coeffs);
}


// ************************************************************************* //
