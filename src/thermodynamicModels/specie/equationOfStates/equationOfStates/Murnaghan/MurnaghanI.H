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
Foam::scalar Foam::Murnaghan<Specie>::p
(
    const scalar& rho,
    const scalar& T
) const
{
    return pRef_ + (pow(rho/rho0_, n_) - 1.0)/(kappa_*n_);
}


template<class Specie>
Foam::scalar Foam::Murnaghan<Specie>::Gamma
(
    const scalar& rho,
    const scalar& T,
    const scalar& cv
) const
{
    return Gamma_ + 1.0;
}


template<class Specie>
Foam::scalar Foam::Murnaghan<Specie>::speedOfSound
(
    const scalar& p,
    const scalar& rho,
    const scalar& T,
    const scalar& cv
) const
{
    return sqrt(pow(rho/rho0_, n_ - 1.0)/(rho0_*n_*kappa_));
}


template<class Specie>
Foam::scalar Foam::Murnaghan<Specie>::dpdv
(
    const scalar& rho,
    const scalar& T
) const
{
    return -rho*pow(rho/rho0_, n_)/kappa_;
}


template<class Specie>
Foam::scalar Foam::Murnaghan<Specie>::dpdT
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::Murnaghan<Specie>::E
(
    const scalar& rho,
    const scalar& T
) const
{
    scalar rhos(max(rho, 1e-10));
    return
       -pRef_/rhos
      + (pow(rho/rho0_, n_) - (n_ - 1.0))/(kappa_*n_*rhos*(n_ - 1.0));
}


template<class Specie>
Foam::scalar Foam::Murnaghan<Specie>::Cv
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::Murnaghan<Specie>::H
(
    const scalar& rho,
    const scalar& T
) const
{
    return pow(rho/rho0_, n_ - 1.0)/(rho0_*kappa_*(n_ - 1.0));
}


template<class Specie>
Foam::scalar Foam::Murnaghan<Specie>::Cp
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::Murnaghan<Specie>::CpMCv
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::Murnaghan<Specie>::S
(
    const scalar& p,
    const scalar& rho,
    const scalar& T
) const
{
    return 0;
}


// ************************************************************************* //
