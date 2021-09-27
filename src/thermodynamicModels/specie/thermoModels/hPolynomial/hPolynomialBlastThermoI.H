/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2021-09-27 Jeff Heylmun:    Modified class for a density based thermodynamic
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

template<class EquationOfState, int PolySize>
inline Foam::hPolynomialThermo<EquationOfState, PolySize>::hPolynomialThermo
(
    const EquationOfState& pt,
    const scalar Hf,
    const scalar Sf,
    const Polynomial<PolySize>& CpCoeffs,
    const typename Polynomial<PolySize>::intPolyType& hCoeffs,
    const Polynomial<PolySize>& sCoeffs
)
:
    EquationOfState(pt),
    Hf_(Hf),
    Sf_(Sf),
    CpCoeffs_(CpCoeffs),
    hCoeffs_(hCoeffs),
    sCoeffs_(sCoeffs)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline Foam::hPolynomialThermo<EquationOfState, PolySize>::hPolynomialThermo
(
    const word& name,
    const hPolynomialThermo& pt
)
:
    EquationOfState(name, pt),
    Hf_(pt.Hf_),
    Sf_(pt.Sf_),
    CpCoeffs_(pt.CpCoeffs_),
    hCoeffs_(pt.hCoeffs_),
    sCoeffs_(pt.sCoeffs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline Foam::scalar
Foam::hPolynomialThermo<EquationOfState, PolySize>::limit
(
    const scalar T,
    const scalar Told
) const
{
    return T;
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<EquationOfState, PolySize>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return CpCoeffs_.value(T) + EquationOfState::Cp(rho, e, T);
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<EquationOfState, PolySize>::Hs
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return hCoeffs_.value(T) + EquationOfState::E(rho, e, T);
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<EquationOfState, PolySize>::Hts
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return hCoeffs_.value(T);
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<EquationOfState, PolySize>::Ha
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Hs(rho, e, T) + Hf();
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<EquationOfState, PolySize>::Hta
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Hts(rho, e, T) + Hf();
}


template<class EquationOfState, int PolySize>
inline Foam::scalar
Foam::hPolynomialThermo<EquationOfState, PolySize>::Hf() const
{
    return Hf_;
}


template<class EquationOfState, int PolySize>
inline Foam::scalar
Foam::hPolynomialThermo<EquationOfState, PolySize>::flameT() const
{
    NotImplemented;
    return 0.0;
}


template<class EquationOfState, int PolySize>
inline Foam::scalar
Foam::hPolynomialThermo<EquationOfState, PolySize>::S
(
    const scalar p,
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return sCoeffs_.value(T) + EquationOfState::S(p, rho, T);
}


template<class EquationOfState, int PolySize>
inline Foam::scalar
Foam::hPolynomialThermo<EquationOfState, PolySize>::dGdT
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class EquationOfState, int PolySize>
inline Foam::scalar
Foam::hPolynomialThermo<EquationOfState, PolySize>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    return CpCoeffs_.derivative(T);
}


template<class EquationOfState, int PolySize>
inline Foam::scalar
Foam::hPolynomialThermo<EquationOfState, PolySize>::dCpdT
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return CpCoeffs_.derivative(T);
}


template<class EquationOfState, int PolySize>
inline Foam::scalar
Foam::hPolynomialThermo<EquationOfState, PolySize>::dCvdT
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    NotImplemented;
    return CpCoeffs_.derivative(T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline void Foam::hPolynomialThermo<EquationOfState, PolySize>::operator+=
(
    const hPolynomialThermo<EquationOfState, PolySize>& pt
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(pt);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        const scalar Y2 = pt.Y()/this->Y();

        Hf_ = Y1*Hf_ + Y2*pt.Hf_;
        Sf_ = Y1*Sf_ + Y2*pt.Sf_;
        CpCoeffs_ = Y1*CpCoeffs_ + Y2*pt.CpCoeffs_;
        hCoeffs_ = Y1*hCoeffs_ + Y2*pt.hCoeffs_;
        sCoeffs_ = Y1*sCoeffs_ + Y2*pt.sCoeffs_;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline Foam::hPolynomialThermo<EquationOfState, PolySize> Foam::operator+
(
    const hPolynomialThermo<EquationOfState, PolySize>& pt1,
    const hPolynomialThermo<EquationOfState, PolySize>& pt2
)
{
    EquationOfState eofs = pt1;
    eofs += pt2;

    if (mag(eofs.Y()) < small)
    {
        return hPolynomialThermo<EquationOfState, PolySize>
        (
            eofs,
            pt1.Hf_,
            pt1.Sf_,
            pt1.CpCoeffs_,
            pt1.hCoeffs_,
            pt1.sCoeffs_
        );
    }
    {
        const scalar Y1 = pt1.Y()/eofs.Y();
        const scalar Y2 = pt2.Y()/eofs.Y();

        return hPolynomialThermo<EquationOfState, PolySize>
        (
            eofs,
            Y1*pt1.Hf_ + Y2*pt2.Hf_,
            Y1*pt1.Sf_ + Y2*pt2.Sf_,
            Y1*pt1.CpCoeffs_ + Y2*pt2.CpCoeffs_,
            Y1*pt1.hCoeffs_ + Y2*pt2.hCoeffs_,
            Y1*pt1.sCoeffs_ + Y2*pt2.sCoeffs_
        );
    }
}


template<class EquationOfState, int PolySize>
inline Foam::hPolynomialThermo<EquationOfState, PolySize> Foam::operator*
(
    const scalar s,
    const hPolynomialThermo<EquationOfState, PolySize>& pt
)
{
    return hPolynomialThermo<EquationOfState, PolySize>
    (
        s*static_cast<const EquationOfState&>(pt),
        pt.Hf_,
        pt.Sf_,
        pt.CpCoeffs_,
        pt.hCoeffs_,
        pt.sCoeffs_
    );
}


template<class EquationOfState, int PolySize>
inline Foam::hPolynomialThermo<EquationOfState, PolySize> Foam::operator==
(
    const hPolynomialThermo<EquationOfState, PolySize>& pt1,
    const hPolynomialThermo<EquationOfState, PolySize>& pt2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(pt1)
     == static_cast<const EquationOfState&>(pt2)
    );

    const scalar Y1 = pt1.Y()/eofs.Y();
    const scalar Y2 = pt2.Y()/eofs.Y();

    return hPolynomialThermo<EquationOfState, PolySize>
    (
        eofs,
        Y2*pt2.Hf_       - Y1*pt1.Hf_,
        Y2*pt2.Sf_       - Y1*pt1.Sf_,
        Y2*pt2.CpCoeffs_ - Y1*pt1.CpCoeffs_,
        Y2*pt2.hCoeffs_  - Y1*pt1.hCoeffs_,
        Y2*pt2.sCoeffs_  - Y1*pt1.sCoeffs_
    );
}


// ************************************************************************* //
