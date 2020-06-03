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

#include "constants.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
Foam::scalar Foam::vanderWaals<Specie>::Pi
(
    const scalar& rho,
    const scalar& e
) const
{
    return
        (1.0 - (gamma_ - 1.0)/stabilise((1.0 - b_*rho), 1e-10))*a_*sqr(rho)
      + ((gamma_ - 1.0)/stabilise((1.0 - b_*rho), 1e-10) + 1.0)*c_;
}


template<class Specie>
Foam::scalar Foam::vanderWaals<Specie>::delta
(
    const scalar& p,
    const scalar& rho,
    const scalar& e
) const
{
    scalar pi
    (
        (gamma_ - 1.0)/stabilise((1.0 - b_*rho), 1e-10)
       *(rho*e + a_*sqr(rho) - c_)
      - (a_*sqr(rho) + c_)
    );
    return
        -b_*(pi + a_*sqr(rho))/(gamma_ - 1.0)
       + ((1.0 - b_*rho)/(gamma_ - 1.0) - 1.0)*2.0*a_*rho;
}


template<class Specie>
Foam::scalar Foam::vanderWaals<Specie>::Gamma
(
    const scalar& rho,
    const scalar& e
) const
{
    return (gamma_ - 1.0)/stabilise((1.0 - b_*rho), 1e-10) + 1.0;
}


template<class Specie>
Foam::scalar Foam::vanderWaals<Specie>::dpdv
(
    const scalar& rho,
    const scalar& e
) const
{
    return
        sqr(rho)/sqr(1.0 - b_*rho)
       *(
            2.0*a_*rho*(b_*rho - 1.0)*(b_*rho + gamma_ - 2.0)
          + b_*c_*(gamma_ - 1.0)
          - a_*b_*sqr(rho)*(gamma_ - 1.0)
          - (gamma_ - 1.0)*e
        );
}


template<class Specie>
Foam::scalar Foam::vanderWaals<Specie>::dpde
(
    const scalar& rho,
    const scalar& e
) const
{
    return (Gamma(rho, e) - 1.0)*rho;
}


template<class Specie>
Foam::scalar Foam::vanderWaals<Specie>::E
(
    const scalar& rho,
    const scalar& e
) const
{
    scalar rhos(max(rho, 1e-10));
    return
      - a_*rho
      + c_*gamma_/rhos
      - (gamma_ - 1.0)/b_
       *(
            a_*log(rhos)
          - (a_ - sqr(b_)*c_)*log(rhos/stabilise((1.0 - b_*rho), 1e-10))
        );
}


template<class Specie>
Foam::scalar Foam::vanderWaals<Specie>::Cv
(
    const scalar& rho,
    const scalar& e
) const
{
    // Check
    return 0;
}


template<class Specie>
Foam::scalar Foam::vanderWaals<Specie>::Cp
(
    const scalar& rho,
    const scalar& e
) const
{
    // Check
    return 0;
}


template<class Specie>
Foam::scalar Foam::vanderWaals<Specie>::H
(
    const scalar& rho,
    const scalar& e
) const
{
    return
        rho/stabilise((1.0 - b_*rho), 1e-10)
       *(
            (gamma_ - 1.0)*e
          + 2.0*a_*b_*rho
          + a_*(gamma_ - 3.0)
          - sqr(b_)*c_*(gamma_ - 1.0)
        )
      - (gamma_ - 1.0)
       *(
            a_*log(rho)
          - (a_ - sqr(b_)*c_)*log(rho/stabilise((1.0 - b_*rho), 1e-10))
        )/b_;
}


template<class Specie>
Foam::scalar Foam::vanderWaals<Specie>::CpMCv
(
    const scalar& rho,
    const scalar& e
) const
{
    return (gamma_ - 1.0)/stabilise((1.0 - b_*rho), 1e-10);
}


template<class Specie>
Foam::scalar Foam::vanderWaals<Specie>::S
(
    const scalar& p,
    const scalar& rho,
    const scalar& T
) const
{
    // Check
    return 0;
}

// ************************************************************************* //
