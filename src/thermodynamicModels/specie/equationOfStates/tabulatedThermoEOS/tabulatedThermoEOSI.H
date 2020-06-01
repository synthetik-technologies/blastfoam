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

template<class Specie>
Foam::scalar Foam::tabulatedThermoEOS<Specie>::p
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return pTable_.lookup(rho, e);
}


template<class Specie>
Foam::scalar Foam::tabulatedThermoEOS<Specie>::dpdv
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return -sqr(rho)*pTable_.dFdX(rho, e);
}


template<class Specie>
Foam::scalar Foam::tabulatedThermoEOS<Specie>::dpdT
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return
        pTable_.dFdX(rho, e)/TTable_.dFdX(rho, e)
      + pTable_.dFdY(rho, e)/TTable_.dFdY(rho, e);
}


template<class Specie>
Foam::scalar Foam::tabulatedThermoEOS<Specie>::dpde
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return pTable_.dFdY(rho, e);
}


template<class Specie>
Foam::scalar Foam::tabulatedThermoEOS<Specie>::Gamma
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return pTable_.dFdY(rho, e)/max(rho, 1e-10);
}


template<class Specie>
Foam::scalar Foam::tabulatedThermoEOS<Specie>::speedOfSound
(
    const scalar& p,
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return
        sqrt
        (
            max
            (
                pTable_.dFdX(rho, e) + Gamma(rho, e, T)*p/max(rho, 1e-10),
                small
            )
        );
}

template<class Specie>
Foam::scalar Foam::tabulatedThermoEOS<Specie>::Es
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return TTable_.reverseLookupY(T, rho);
}


template<class Specie>
Foam::scalar Foam::tabulatedThermoEOS<Specie>::Ea
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return TTable_.reverseLookupY(T, rho);
}

template<class Specie>
Foam::scalar Foam::tabulatedThermoEOS<Specie>::Cv
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return 1.0/max(TTable_.dFdY(rho, e), small);
}


template<class Specie>
Foam::scalar Foam::tabulatedThermoEOS<Specie>::Cp
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    scalar cv(Cv(rho, e, T));
    scalar cSqr
    (
        sqr
        (
            speedOfSound
            (
                p(rho, e, T),
                rho,
                e,
                T
            )
        )
    );
    return cv/(1.0 - sqr(Gamma(rho, e, T))*cv*T/cSqr);
}


template<class Specie>
Foam::scalar Foam::tabulatedThermoEOS<Specie>::CpMCv
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return Cp(rho, e, T) - Cv(rho, e, T);
}

template<class Specied>
Foam::scalar Foam::tabulatedThermoEOS<Specied>::CpByCv
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return Cp(rho, e, T)/Cv(rho, e, T);
}


template<class Specied>
Foam::scalar Foam::tabulatedThermoEOS<Specied>::TRhoE
(
    const scalar& T0,
    const scalar& rho,
    const scalar& e
) const
{
    return TTable_.lookup(rho, e);
}


template<class Specied>
Foam::scalar Foam::tabulatedThermoEOS<Specied>::initializeEnergy
(
    const scalar& p,
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return pTable_.reverseLookupY(p, rho);
}

template<class Specied>
Foam::scalar Foam::tabulatedThermoEOS<Specied>::initializeRho
(
    const scalar& p,
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return rho;
}
// ************************************************************************* //
