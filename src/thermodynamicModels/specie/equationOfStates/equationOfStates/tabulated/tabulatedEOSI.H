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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::autoPtr<Foam::tabulatedEOS<Specie>>
Foam::tabulatedEOS<Specie>::clone() const
{
    return autoPtr<tabulatedEOS<Specie>>
    (
        new tabulatedEOS<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::tabulatedEOS<Specie>>
Foam::tabulatedEOS<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<tabulatedEOS<Specie>>
    (
        new tabulatedEOS<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
Foam::scalar Foam::tabulatedEOS<Specie>::p
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const bool limit
) const
{
    return limit ? max(pRhoT(rho, T), 0.0) : pRhoT(rho, T);
}


template<class Specie>
Foam::scalar Foam::tabulatedEOS<Specie>::pRhoT
(
    const scalar rho,
    const scalar T
) const
{
    return pTable_.lookup(rho, T);
}


template<class Specie>
Foam::scalar Foam::tabulatedEOS<Specie>::Gamma
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar cv
) const
{
    return max(dpdT(rho, e, T)/(cv*max(rho, 1e-10)), small) + 1.0;
}


template<class Specie>
Foam::scalar Foam::tabulatedEOS<Specie>::cSqr
(
    const scalar p,
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar cv
) const
{
    return
        (sqr(dpdT(rho, e, T))*T/cv - dpdv(rho, e, T))
       /sqr(max(rho, 1e-10));
}


template<class Specie>
Foam::scalar Foam::tabulatedEOS<Specie>::dpdv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return -sqr(rho)*pTable_.dFdX(rho, e);
}


template<class Specie>
Foam::scalar Foam::tabulatedEOS<Specie>::dpde
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    NotImplemented;
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::tabulatedEOS<Specie>::dpdT
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return pTable_.dFdY(rho, T);
}


template<class Specie>
Foam::scalar Foam::tabulatedEOS<Specie>::E
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class Specie>
Foam::scalar Foam::tabulatedEOS<Specie>::Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class Specie>
Foam::scalar Foam::tabulatedEOS<Specie>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class Specie>
Foam::scalar Foam::tabulatedEOS<Specie>::H
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class Specie>
Foam::scalar Foam::tabulatedEOS<Specie>::CpMCv
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar CpCv
) const
{
    return pTable_.dFdY(rho, T)/max(rho, 1e-10);
}


template<class Specie>
Foam::scalar Foam::tabulatedEOS<Specie>::S
(
    const scalar p,
    const scalar rho,
    const scalar T
) const
{
    NotImplemented;
    return 0.0;
}


// ************************************************************************* //
