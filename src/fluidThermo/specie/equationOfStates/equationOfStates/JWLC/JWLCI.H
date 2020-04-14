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

#include "JWLC.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
Foam::scalar Foam::JWLC<Specie>::Gamma
(
    const scalar& rho,
    const scalar& T,
    const scalar& cv
) const
{
    return omega_ + 1.0;
}


template<class Specie>
Foam::scalar Foam::JWLC<Specie>::p
(
    const scalar& rho,
    const scalar& T
) const
{
    scalar V(rho0_/max(rho, 1e-10));
    return A_*exp(-R1_*V) + B_*exp(-R2_*V) + C_/pow(V, omega_ + 1.0);
}


template<class Specie>
Foam::scalar Foam::JWLC<Specie>::speedOfSound
(
    const scalar& p,
    const scalar& rho,
    const scalar& T,
    const scalar& cv
) const
{
    scalar V(rho0_/max(rho, 1e-10));
    return
        sqr(V)/rho0_
       *(
            A_*R1_*exp(-R1_*V)
          + B_*R2_*exp(-R2_*V)
          + C_*(1.0 + omega_)/pow(V, 2.0 + omega_)
        );
}


template<class Specie>
Foam::scalar Foam::JWLC<Specie>::dpdv
(
    const scalar& rho,
    const scalar& T
) const
{
    scalar V(rho0_/max(rho, 1e-10));
    return
        rho0_
       *(
            A_*R1_*exp(-R1_*V)
          + B_*R2_*exp(-R2_*V)
          + C_*(1.0 + omega_)/pow(V, 2.0 + omega_)
        );
}


template<class Specie>
Foam::scalar Foam::JWLC<Specie>::dpdT
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::JWLC<Specie>::E
(
    const scalar& rho,
    const scalar& T
) const
{
    scalar V(rho0_/max(rho, 1e-10));
    return
        (
            A_/R1_*exp(-R1_*V)
          + B_/R2_*exp(-R2_*V)
          + C_/(omega_*pow(V, omega_))
        )/rho0_;
}


template<class Specie>
Foam::scalar Foam::JWLC<Specie>::Cv
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0;
}


template<class Specie>
Foam::scalar Foam::JWLC<Specie>::Cp
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::JWLC<Specie>::H
(
    const scalar& rho,
    const scalar& T
) const
{
    scalar rhos(max(rho, 1e-10));
    scalar V(rho0_/rhos);
    return
      - (
            A_*R1_*rho0_*exp(-R1_*V)
          + B_*R2_*rho0_*exp(-R2_*V)
        )/(2.0*sqr(rhos))
      - C_*(omega_ + 1.0)*pow(V, -omega_)/(2.0*rho0_);
}


template<class Specie>
Foam::scalar Foam::JWLC<Specie>::CpMCv
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::JWLC<Specie>::S
(
    const scalar& p,
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}

// ************************************************************************* //
