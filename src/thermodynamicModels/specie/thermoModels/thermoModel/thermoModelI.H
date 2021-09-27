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
#include "IOmanip.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
inline Foam::thermoModel<ThermoType>::thermoModel
(
    const ThermoType& sp
)
:
    ThermoType(sp),
    relTol_(1e-4),
    absTol_(1e-2),
    maxIter_(100)
{
    this->eBased_ = !ThermoType::enthalpy();
}


template<class ThermoType>
inline Foam::thermoModel<ThermoType>::thermoModel
(
    const thermoType& thermo
)
:
    ThermoType(thermo),
    relTol_(thermo.relTol_),
    absTol_(thermo.absTol_),
    maxIter_(thermo.maxIter_)
{
    this->eBased_ = !ThermoType::enthalpy();
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::T
(
    const scalar f,
    const scalar e,
    const scalar rho,
    const scalar T0,
    scalar (thermoModel<ThermoType>::*F)(const scalar, const scalar, const scalar) const,
    scalar (thermoModel<ThermoType>::*dFdT)(const scalar, const scalar, const scalar)
        const,
    scalar (thermoModel<ThermoType>::*limit)(const scalar, const scalar) const,
    const bool diagnostics
) const
{
    if (rho < small)
    {
        return T0;
    }

    scalar Test = (this->*limit)(T0, T0);
    scalar Tnew = Test;

    int    iter = 0;
    scalar relError = great;
    scalar absError = great;

    if (diagnostics)
    {
        const unsigned int width = IOstream::defaultPrecision() + 8;

        InfoInFunction
            << "Energy -> temperature conversion failed to converge:" << endl;
        Pout<< setw(width) << "iter"
            << setw(width) << "Test"
            << setw(width) << "e/h Actual"
            << setw(width) << "e/h Test"
            << setw(width) << "Cv/Cp"
            << setw(width) << "Tnew"
            << endl;
    }

    do
    {
        Test = Tnew;
        Tnew =
            (this->*limit)
            (
                Test
              - ((this->*F)(rho, f, Test) - f)/(this->*dFdT)(rho, f, Test),
                Test
            );

        absError = mag((Tnew - Test));
        relError = absError/max(Tnew, small);

        if (diagnostics)
        {
            const unsigned int width = IOstream::defaultPrecision() + 8;

            Pout<< setw(width) << iter
                << setw(width) << Test
                << setw(width) << f
                << setw(width) << ((this->*F)(rho, f, Test))
                << setw(width) << ((this->*dFdT)(rho, f, Test))
                << setw(width) << Tnew
                << endl;
        }

        if (iter++ > maxIter_)
        {
            if (!diagnostics)
            {
                T(f, e, rho, T0, F, dFdT, limit, true);
            }

            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << maxIter_
                << abort(FatalError);
        }
    } while (relError > relTol_ && absError > absTol_);

    return Tnew;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
inline Foam::thermoModel<ThermoType>::thermoModel
(
    const word& name,
    const ThermoType& sp
)
:
    ThermoType(name, sp),
    relTol_(1e-2),
    absTol_(1e-4),
    maxIter_(100)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
Foam::scalar Foam::thermoModel<ThermoType>::Gamma
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return ThermoType::Gamma(rho, e, T, ThermoType::Cv(rho, e, T));
}


template<class ThermoType>
Foam::scalar Foam::thermoModel<ThermoType>::p
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return ThermoType::p(rho, e, T, true);
}


template<class ThermoType>
Foam::scalar Foam::thermoModel<ThermoType>::dpdRho
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return -ThermoType::dpdv(rho, e, T)/sqr(max(rho, 1e-10));
}


template<class ThermoType>
Foam::scalar Foam::thermoModel<ThermoType>::dpde
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (ThermoType::temperatureBased())
    {
        return ThermoType::dpdT(rho, e, T)/ThermoType::Cv(rho, e, T);
    }
    return ThermoType::dpde(rho, e, T);
}


template<class ThermoType>
Foam::scalar Foam::thermoModel<ThermoType>::cSqr
(
    const scalar p,
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return ThermoType::cSqr(p, rho, e, T, ThermoType::Cv(rho, e, T));
}


template<class ThermoType>
Foam::scalar Foam::thermoModel<ThermoType>::CpByCv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return ThermoType::Cp(rho, e, T)/ThermoType::Cv(rho, e, T);
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::G
(
    const scalar p,
    const scalar rho,
    const scalar T
) const
{
    scalar e = this->Es(rho, 1000.0, T);
    return this->Ha(rho, e, T) - T*this->S(p, rho, e, T);
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::A
(
    const scalar p,
    const scalar rho,
    const scalar T
) const
{
    scalar e = this->Es(rho, 1000.0, T);
    return this->Ea(rho, e, T) - T*this->S(p, rho, e, T);
}


template<class ThermoType>
inline Foam::scalar
Foam::thermoModel<ThermoType>::K
(
    const scalar p,
    const scalar rho,
    const scalar T
) const
{
    scalar arg = -this->Y()*this->G(Pstd, rho, T)/(RR*T);

    if (arg < 600)
    {
        return exp(arg);
    }
    else
    {
        return rootVGreat;
    }
}


template<class ThermoType>
inline Foam::scalar
Foam::thermoModel<ThermoType>::Kp
(
    const scalar p,
    const scalar rho,
    const scalar T
) const
{
    return K(p, rho, T);
}


template<class ThermoType>
inline Foam::scalar
Foam::thermoModel<ThermoType>::Kc
(
    const scalar p,
    const scalar rho,
    const scalar T
) const
{
    const scalar nm = this->Y()/this->W();

    if (equal(nm, small))
    {
        return Kp(p, rho, T);
    }
    else
    {
        return Kp(p, rho, T)*pow(Pstd/(RR*T), nm);
    }
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::Kx
(
    const scalar p,
    const scalar rho,
    const scalar T
) const
{
    const scalar nm = this->Y()/this->W();

    if (equal(nm, small))
    {
        return Kp(p, rho, T);
    }
    else
    {
        return Kp(p, rho, T)*pow(Pstd/p, nm);
    }
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::Kn
(
    const scalar p,
    const scalar rho,
    const scalar T,
    const scalar n
) const
{
    const scalar nm = this->Y()/this->W();

    if (equal(nm, small))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, rho, T)*pow(n*Pstd/p, nm);
    }
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return this->Cp(rho, e, T)*this->W();
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::ha
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return this->Ha(rho, e, T)*this->W();
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::hs
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return this->Hs(rho, e, T)*this->W();
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::hc() const
{
    return this->Hf()*this->W();
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::s
(
    const scalar p,
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return this->S(p, rho, e, T)*this->W();
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::he
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return this->es(rho, e, T);
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return this->Cv(rho, e, T)*this->W();
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::es
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return this->Es(rho, e, T)*this->W();
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::ea
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return this->Ea(rho, e, T)*this->W();
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::g
(
    const scalar p,
    const scalar rho,
    const scalar T
) const
{
    return this->G(p, rho, T)*this->W();
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::a
(
    const scalar p,
    const scalar rho,
    const scalar T
) const
{
    return this->A(p, rho, T)*this->W();
}



template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::dKcdTbyKc
(
    const scalar p,
    const scalar T
) const
{
    const scalar nm = this->Y()/this->W();

    if (equal(nm, small))
    {
        return -this->dGdT(Pstd, T)*this->Y()/RR;
    }
    else
    {
        return -(nm/T + this->dGdT(Pstd, T)*this->Y()/RR);
    }
}

template<class ThermoType>
inline Foam::scalar
Foam::thermoModel<ThermoType>::dcpdT(const scalar p, const scalar T) const
{
    return this->dCpdT(p, T)*this->W();
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::initializeEnergy
(
    const scalar p,
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (rho < small)
    {
        return 0.0;
    }

    scalar Eest = ThermoType::Es(rho, e, T);

    if (ThermoType::temperatureBased())
    {
        return Eest;
    }

    if (mag(ThermoType::dpde(rho, e, T)) < small)
    {
        return Eest;
    }

    scalar Enew = Eest;
    scalar pNew;
    int    iter = 0;
    do
    {
        Eest = Enew;
        pNew = ThermoType::p(rho, Eest, T, false); // Do not limit p
        Enew -=
            (pNew - p)/stabilise(ThermoType::dpde(rho, Eest, T), small);

        if (iter++ > 100)
        {
            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << 100
                << abort(FatalError);
        }

    } while ((mag(pNew - p)/p > relTol_) && (iter++ < maxIter_));
    return Enew;
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::rhoPT
(
    const scalar p,
    const scalar T
) const
{
    //- Simple method to calculate initial density
    //  Should be modified to solve the 2D problem for
    //  density and internal energy
    scalar Rhoest = 1.0;
    scalar Rhonew = Rhoest;
    scalar pNew = p;
    scalar E = ThermoType::Es(Rhoest, 0.0, T); //- Initial guess

    int    iter = 0;
    do
    {
        Rhoest = Rhonew;

        scalar dpdRho(-ThermoType::dpdv(Rhoest, E, T)/sqr(max(Rhoest, 1e-10)));
        pNew = ThermoType::p(Rhoest, E, T, false); // Do not limit p
        Rhonew =
            Rhoest - (pNew - p)/stabilise(dpdRho, small);
        E = ThermoType::Es(Rhoest, E, T);

        if (Rhonew < 1e-10)
        {
            Rhonew = Rhoest/2.0;
        }
    } while ((mag(pNew - p)/p > relTol_) && (iter++ < maxIter_));

    return Rhonew;
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::TRhoE
(
    const scalar T0,
    const scalar rho,
    const scalar e
) const
{
    return T
    (
        e,
        e,
        rho,
        T0,
        &thermoModel<ThermoType>::Es,
        &thermoModel<ThermoType>::Cv,
        &thermoModel<ThermoType>::limit
    );
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::THs
(
    const scalar hs,
    const scalar p,
    const scalar rho,
    const scalar T0
) const
{
    return T
    (
        hs,
        hs - p/rho,
        rho,
        T0,
        &thermoModel<ThermoType>::Hs,
        &thermoModel<ThermoType>::Cp,
        &thermoModel<ThermoType>::limit
    );
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::THa
(
    const scalar ha,
    const scalar p,
    const scalar rho,
    const scalar T0
) const
{
    return T
    (
        ha,
        ha - p/rho - this->Hf(),
        rho,
        T0,
        &thermoModel<ThermoType>::Ha,
        &thermoModel<ThermoType>::Cp,
        &thermoModel<ThermoType>::limit
    );
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::TEs
(
    const scalar es,
    const scalar p,
    const scalar rho,
    const scalar T0
) const
{
    return T
    (
        es,
        es,
        rho,
        T0,
        &thermoModel<ThermoType>::Es,
        &thermoModel<ThermoType>::Cv,
        &thermoModel<ThermoType>::limit
    );
}


template<class ThermoType>
inline Foam::scalar Foam::thermoModel<ThermoType>::TEa
(
    const scalar ea,
    const scalar p,
    const scalar rho,
    const scalar T0
) const
{
    return T
    (
        ea,
        ea - this->Hf(),
        rho,
        T0,
        &thermoModel<ThermoType>::Ea,
        &thermoModel<ThermoType>::Cv,
        &thermoModel<ThermoType>::limit
    );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class ThermoType>
inline void Foam::thermoModel<ThermoType>::operator=
(
    const thermoModel<ThermoType>& st
)
{
    ThermoType::operator=(st);
}


template<class ThermoType>
inline void Foam::thermoModel<ThermoType>::operator+=
(
    const thermoModel<ThermoType>& st
)
{
    ThermoType::operator+=(st);
}


template<class ThermoType>
inline void Foam::thermoModel<ThermoType>::operator*=(const scalar s)
{
    ThermoType::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class ThermoType>
inline Foam::thermoModel<ThermoType> Foam::operator+
(
    const thermoModel<ThermoType>& st1,
    const thermoModel<ThermoType>& st2
)
{
    return thermoModel<ThermoType>
    (
        static_cast<const ThermoType&>(st1)
      + static_cast<const ThermoType&>(st2)
    );
}


template<class ThermoType>
inline Foam::thermoModel<ThermoType> Foam::operator*
(
    const scalar s,
    const thermoModel<ThermoType>& st
)
{
    return thermoModel<ThermoType>
    (
        s*static_cast<const ThermoType&>(st)
    );
}


template<class ThermoType>
inline Foam::thermoModel<ThermoType> Foam::operator==
(
    const thermoModel<ThermoType>& st1,
    const thermoModel<ThermoType>& st2
)
{
    return thermoModel<ThermoType>
    (
        static_cast<const ThermoType&>(st1)
     == static_cast<const ThermoType&>(st2)
    );
}


// ************************************************************************* //
