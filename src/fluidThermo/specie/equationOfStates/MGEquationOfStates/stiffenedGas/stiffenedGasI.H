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

#include "stiffenedGas.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
Foam::scalar Foam::stiffenedGas<Specie>::Pi
(
    const scalar& rho,
    const scalar& e
) const
{
    return gamma_*a_;
}


template<class Specie>
Foam::scalar Foam::stiffenedGas<Specie>::delta
(
    const scalar& p,
    const scalar& rho,
    const scalar& e
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::stiffenedGas<Specie>::Gamma
(
    const scalar& rho,
    const scalar& e
) const
{
    return gamma_;
}


template<class Specie>
Foam::scalar Foam::stiffenedGas<Specie>::dpdv
(
    const scalar& rho,
    const scalar& e
) const
{
    return -sqr(rho)*(gamma_ - 1.0)*e;
}


template<class Specie>
Foam::scalar Foam::stiffenedGas<Specie>::dpde
(
    const scalar& rho,
    const scalar& e
) const
{
    return (gamma_ - 1.0)*rho;
}


template<class Specie>
Foam::scalar Foam::stiffenedGas<Specie>::E
(
    const scalar& rho,
    const scalar& e
) const
{
    return 0.0;//gamma_*a_/max(rho, 1.0e-10);
}


template<class Specie>
Foam::scalar Foam::stiffenedGas<Specie>::Cv
(
    const scalar& rho,
    const scalar& e
) const
{
    // Check
    return 0;
}


template<class Specie>
Foam::scalar Foam::stiffenedGas<Specie>::Cp
(
    const scalar& rho,
    const scalar& e
) const
{
    // Check
    return 0;
}


template<class Specie>
Foam::scalar Foam::stiffenedGas<Specie>::H
(
    const scalar& rho,
    const scalar& e
) const
{
    return 0;
}


template<class Specie>
Foam::scalar Foam::stiffenedGas<Specie>::CpMCv
(
    const scalar& rho,
    const scalar& e
) const
{
    return this->R();
}


template<class Specie>
Foam::scalar Foam::stiffenedGas<Specie>::S
(
    const scalar& p,
    const scalar& rho,
    const scalar& T
) const
{
    scalar rho0((Pstd - gamma_*a_)/this->R()*Tstd);
    return this->R()*(log(rho0) - log(rho));
}

// ************************************************************************* //
