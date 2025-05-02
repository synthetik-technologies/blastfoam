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
inline Foam::janaf9Thermo<EquationOfState>::janaf9Thermo
(
    const EquationOfState& st,
    const scalar Tlow,
    const scalar Thigh,
    const scalar Tcommon,
    const typename janaf9Thermo<EquationOfState>::coeffArray& highCpCoeffs,
    const typename janaf9Thermo<EquationOfState>::coeffArray& lowCpCoeffs,
    const bool convertCoeffs
)
:
    EquationOfState(st),
    Tlow_(Tlow),
    Thigh_(Thigh),
    Tcommon_(Tcommon)
{
    if (convertCoeffs)
    {
        for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
        {
            highCpCoeffs_[coefLabel] = highCpCoeffs[coefLabel]*this->R();
            lowCpCoeffs_[coefLabel] = lowCpCoeffs[coefLabel]*this->R();
        }
    }
    else
    {
        for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
        {
            highCpCoeffs_[coefLabel] = highCpCoeffs[coefLabel];
            lowCpCoeffs_[coefLabel] = lowCpCoeffs[coefLabel];
        }
    }
}


template<class EquationOfState>
inline const typename Foam::janaf9Thermo<EquationOfState>::coeffArray&
Foam::janaf9Thermo<EquationOfState>::coeffs
(
    const scalar T
) const
{
    if (T < Tcommon_)
    {
        return lowCpCoeffs_;
    }
    else
    {
        return highCpCoeffs_;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::janaf9Thermo<EquationOfState>::janaf9Thermo
(
    const word& name,
    const janaf9Thermo& jt
)
:
    EquationOfState(name, jt),
    Tlow_(jt.Tlow_),
    Thigh_(jt.Thigh_),
    Tcommon_(jt.Tcommon_)
{
    for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
    {
        highCpCoeffs_[coefLabel] = jt.highCpCoeffs_[coefLabel];
        lowCpCoeffs_[coefLabel] = jt.lowCpCoeffs_[coefLabel];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::janaf9Thermo<EquationOfState>::limit
(
    const scalar T,
    const scalar Told
) const
{
    if (T < Tlow_ || T > Thigh_)
    {
        if (janaf9Thermo<EquationOfState>::debug)
        {
            WarningInFunction
                << "attempt to use janaf9Thermo<EquationOfState>"
                << " out of temperature range "
                << Tlow_ << " -> " << Thigh_ << ";  T = " << T
                << endl;
        }

        return min(max(T, Tlow_), Thigh_);
    }
    else
    {
        return T;
    }
}


template<class EquationOfState>
inline Foam::scalar Foam::janaf9Thermo<EquationOfState>::Tlow() const
{
    return Tlow_;
}


template<class EquationOfState>
inline Foam::scalar Foam::janaf9Thermo<EquationOfState>::Thigh() const
{
    return Thigh_;
}


template<class EquationOfState>
inline Foam::scalar Foam::janaf9Thermo<EquationOfState>::Tcommon() const
{
    return Tcommon_;
}


template<class EquationOfState>
inline const typename Foam::janaf9Thermo<EquationOfState>::coeffArray&
Foam::janaf9Thermo<EquationOfState>::highCpCoeffs() const
{
    return highCpCoeffs_;
}


template<class EquationOfState>
inline const typename Foam::janaf9Thermo<EquationOfState>::coeffArray&
Foam::janaf9Thermo<EquationOfState>::lowCpCoeffs() const
{
    return lowCpCoeffs_;
}


template<class EquationOfState>
inline Foam::scalar Foam::janaf9Thermo<EquationOfState>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Cpt(T) + EquationOfState::Cp(rho, e, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::janaf9Thermo<EquationOfState>::Ha
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Hat(T) + EquationOfState::H(rho, e, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::janaf9Thermo<EquationOfState>::Hs
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Ha(rho, e, T) - Hf();
}

template<class EquationOfState>
inline Foam::scalar Foam::janaf9Thermo<EquationOfState>::Cpt
(
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    const scalar rT = 1.0/T;
    return
        rT*(rT*a[0] + a[1])
      + a[2]+ T*(a[3] + T*(a[4] + T*(a[5] + T*a[6])));
}

template<class EquationOfState>
inline Foam::scalar Foam::janaf9Thermo<EquationOfState>::Hat
(
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    return
      - a[0]/T + a[1]*log(T)
      + a[7]
      + T*(a[2] + T*(0.5*a[3] + T*(a[4]/3.0 + T*(0.25*a[5] + T*0.2*a[6]))));
}

template<class EquationOfState>
inline Foam::scalar Foam::janaf9Thermo<EquationOfState>::Hst
(
    const scalar T
) const
{
    return Hat(T) - Hf();
}


template<class EquationOfState>
inline Foam::scalar Foam::janaf9Thermo<EquationOfState>::Hf() const
{
    return Hat(Tstd);

}


template<class EquationOfState>
inline Foam::scalar Foam::janaf9Thermo<EquationOfState>::S
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    const scalar rT = 1.0/T;
    return
    (
      - rT*(0.5*a[0]*rT + a[1])
      + a[2]*log(T) + T*(a[3] + T*(0.5*a[4] + T*(a[5]/3.0 + T*0.25*a[6])))
      + a[8]
    ) + EquationOfState::S(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::janaf9Thermo<EquationOfState>::Gstd
(
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    const scalar lnT = log(T);
    return
    (
      - 0.5*a[0]/T
      + a[1]*(1.0 + lnT)
      + (
            a[2]*(1.0 - lnT)
          - 0.5*(a[3] + T*(a[4]/3.0 + T*(a[5]/6.0 + T*a[6]/10.0)))
          - a[8]
        )*T
      + a[7]
    );
}


template<class EquationOfState>
inline Foam::scalar Foam::janaf9Thermo<EquationOfState>::dCpdT
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    const scalar rT = 1.0/T;
    return
      - (2.0*a[0]*rT + a[1])*sqr(rT)
      + a[3] + T*(2.0*a[4] + T*(3.0*a[5] + 4.0*T*a[6]));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::janaf9Thermo<EquationOfState>::operator+=
(
    const janaf9Thermo<EquationOfState>& jt
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(jt);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        const scalar Y2 = jt.Y()/this->Y();

        Tlow_ = max(Tlow_, jt.Tlow_);
        Thigh_ = min(Thigh_, jt.Thigh_);

        if
        (
            janaf9Thermo<EquationOfState>::debug
         && notEqual(Tcommon_, jt.Tcommon_)
        )
        {
            FatalErrorInFunction
                << "Tcommon " << Tcommon_ << " for "
                << (this->name().size() ? this->name() : "others")
                << " != " << jt.Tcommon_ << " for "
                << (jt.name().size() ? jt.name() : "others")
                << exit(FatalError);
        }

        for
        (
            label coefLabel=0;
            coefLabel<janaf9Thermo<EquationOfState>::nCoeffs_;
            coefLabel++
        )
        {
            highCpCoeffs_[coefLabel] =
                Y1*highCpCoeffs_[coefLabel]
              + Y2*jt.highCpCoeffs_[coefLabel];

            lowCpCoeffs_[coefLabel] =
                Y1*lowCpCoeffs_[coefLabel]
              + Y2*jt.lowCpCoeffs_[coefLabel];
        }
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::janaf9Thermo<EquationOfState> Foam::operator+
(
    const janaf9Thermo<EquationOfState>& jt1,
    const janaf9Thermo<EquationOfState>& jt2
)
{
    EquationOfState eofs = jt1;
    eofs += jt2;

    if (mag(eofs.Y()) < small)
    {
        return janaf9Thermo<EquationOfState>
        (
            eofs,
            jt1.Tlow_,
            jt1.Thigh_,
            jt1.Tcommon_,
            jt1.highCpCoeffs_,
            jt1.lowCpCoeffs_
        );

    }
    else
    {
        const scalar Y1 = jt1.Y()/eofs.Y();
        const scalar Y2 = jt2.Y()/eofs.Y();

        typename janaf9Thermo<EquationOfState>::coeffArray highCpCoeffs;
        typename janaf9Thermo<EquationOfState>::coeffArray lowCpCoeffs;

        for
        (
            label coefLabel=0;
            coefLabel<janaf9Thermo<EquationOfState>::nCoeffs_;
            coefLabel++
        )
        {
            highCpCoeffs[coefLabel] =
                Y1*jt1.highCpCoeffs_[coefLabel]
              + Y2*jt2.highCpCoeffs_[coefLabel];

            lowCpCoeffs[coefLabel] =
                Y1*jt1.lowCpCoeffs_[coefLabel]
              + Y2*jt2.lowCpCoeffs_[coefLabel];
        }

        if
        (
            janaf9Thermo<EquationOfState>::debug
         && notEqual(jt1.Tcommon_, jt2.Tcommon_)
        )
        {
            FatalErrorInFunction
                << "Tcommon " << jt1.Tcommon_ << " for "
                << (jt1.name().size() ? jt1.name() : "others")
                << " != " << jt2.Tcommon_ << " for "
                << (jt2.name().size() ? jt2.name() : "others")
                << exit(FatalError);
        }

        return janaf9Thermo<EquationOfState>
        (
            eofs,
            max(jt1.Tlow_, jt2.Tlow_),
            min(jt1.Thigh_, jt2.Thigh_),
            jt1.Tcommon_,
            highCpCoeffs,
            lowCpCoeffs
        );
    }
}


template<class EquationOfState>
inline Foam::janaf9Thermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const janaf9Thermo<EquationOfState>& jt
)
{
    return janaf9Thermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(jt),
        jt.Tlow_,
        jt.Thigh_,
        jt.Tcommon_,
        jt.highCpCoeffs_,
        jt.lowCpCoeffs_
    );
}


template<class EquationOfState>
inline Foam::janaf9Thermo<EquationOfState> Foam::operator==
(
    const janaf9Thermo<EquationOfState>& jt1,
    const janaf9Thermo<EquationOfState>& jt2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(jt1)
     == static_cast<const EquationOfState&>(jt2)
    );

    const scalar Y1 = jt2.Y()/eofs.Y();
    const scalar Y2 = jt1.Y()/eofs.Y();

    typename janaf9Thermo<EquationOfState>::coeffArray highCpCoeffs;
    typename janaf9Thermo<EquationOfState>::coeffArray lowCpCoeffs;

    for
    (
        label coefLabel=0;
        coefLabel<janaf9Thermo<EquationOfState>::nCoeffs_;
        coefLabel++
    )
    {
        highCpCoeffs[coefLabel] =
            Y1*jt2.highCpCoeffs_[coefLabel]
          - Y2*jt1.highCpCoeffs_[coefLabel];

        lowCpCoeffs[coefLabel] =
            Y1*jt2.lowCpCoeffs_[coefLabel]
          - Y2*jt1.lowCpCoeffs_[coefLabel];
    }

    if
    (
        janaf9Thermo<EquationOfState>::debug
     && notEqual(jt2.Tcommon_, jt1.Tcommon_)
    )
    {
        FatalErrorInFunction
            << "Tcommon " << jt2.Tcommon_ << " for "
            << (jt2.name().size() ? jt2.name() : "others")
            << " != " << jt1.Tcommon_ << " for "
            << (jt1.name().size() ? jt1.name() : "others")
            << exit(FatalError);
    }

    return janaf9Thermo<EquationOfState>
    (
        eofs,
        max(jt2.Tlow_, jt1.Tlow_),
        min(jt2.Thigh_, jt1.Thigh_),
        jt2.Tcommon_,
        highCpCoeffs,
        lowCpCoeffs
    );
}

// ************************************************************************* //
