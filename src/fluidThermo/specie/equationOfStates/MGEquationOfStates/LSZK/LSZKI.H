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
Foam::scalar Foam::LSZK<Specie>::Pi
(
    const scalar& rho,
    const scalar& e
) const
{
    return -a_*pow(rho, b_);
}


template<class Specie>
Foam::scalar Foam::LSZK<Specie>::delta
(
    const scalar& p,
    const scalar& rho,
    const scalar& e
) const
{
    return -a_*b_*pow(rho, b_ - 1.0)/(gamma_ - 1.0);
}


template<class Specie>
Foam::scalar Foam::LSZK<Specie>::Gamma
(
    const scalar& rho,
    const scalar& e
) const
{
    return gamma_;
}


template<class Specie>
Foam::scalar Foam::LSZK<Specie>::dpdv
(
    const scalar& rho,
    const scalar& e
) const
{
    return -sqr(rho)*(gamma_ - 1.0)*e - a_*b_*pow(rho, b_ + 1.0);
}


template<class Specie>
Foam::scalar Foam::LSZK<Specie>::dpde
(
    const scalar& rho,
    const scalar& e
) const
{
    return (gamma_ - 1.0)*rho;
}


template<class Specie>
Foam::scalar Foam::LSZK<Specie>::E
(
    const scalar& rho,
    const scalar& e
) const
{
    return a_*pow(rho, b_ - 1.0)/(b_ - 1.0);
}


template<class Specie>
Foam::scalar Foam::LSZK<Specie>::Cv
(
    const scalar& rho,
    const scalar& e
) const
{
    // Check
    return 0;
}


template<class Specie>
Foam::scalar Foam::LSZK<Specie>::Cp
(
    const scalar& rho,
    const scalar& e
) const
{
    // Check
    return 0;
}


template<class Specie>
Foam::scalar Foam::LSZK<Specie>::H
(
    const scalar& rho,
    const scalar& e
) const
{
    return a_*log(rho);
}


template<class Specie>
Foam::scalar Foam::LSZK<Specie>::CpMCv
(
    const scalar& rho,
    const scalar& e
) const
{
    return this->R();
}


template<class Specie>
Foam::scalar Foam::LSZK<Specie>::S
(
    const scalar& p,
    const scalar& rho,
    const scalar& T
) const
{
    // check
    return -this->R()*log(p/Foam::constant::thermodynamic::Pstd);
}


// ************************************************************************* //
