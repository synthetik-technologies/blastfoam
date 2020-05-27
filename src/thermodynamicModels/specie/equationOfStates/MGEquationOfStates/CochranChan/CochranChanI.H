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
Foam::scalar Foam::CochranChan<Specie>::Gamma
(
    const scalar& rho,
    const scalar& e
) const
{
    return Gamma0_ + 1.0;
}

template<class Specie>
Foam::scalar Foam::CochranChan<Specie>::Pi
(
    const scalar& rho,
    const scalar& e
) const
{
    scalar rhos(max(rho, 1e-10));
    return
        Gamma0_*rho
       *(
          - A_/((1.0 - Epsilon1_)*rho0_)*(pow(rho0_/rhos, 1.0 - Epsilon1_) - 1.0)
          + B_/((1.0 - Epsilon2_)*rho0_)*(pow(rho0_/rhos, 1.0 - Epsilon2_) - 1.0)
        )
      - A_*pow(rho0_/rhos, -Epsilon1_) + B_*pow(rho0_/rhos, -Epsilon2_);
}


template<class Specie>
Foam::scalar Foam::CochranChan<Specie>::delta
(
    const scalar& p,
    const scalar& rho,
    const scalar& e
) const
{
    scalar rhos(max(rho, 1e-10));
    return
      (
           -A_
           *(
               Epsilon1_*pow(rho0_/rhos, -Epsilon1_)
              *(Epsilon1_ - Gamma0_ - 1.0)/rhos
             + Gamma0_/rho0_
            )/(Epsilon1_ - 1.0)
          + B_
           *(
               Epsilon2_*pow(rho0_/rhos, -Epsilon2_)
              *(Epsilon2_ - Gamma0_ - 1.0)/rhos
             + Gamma0_/rho0_
            )/(Epsilon2_ - 1.0)
        )/Gamma0_;
}


template<class Specie>
Foam::scalar Foam::CochranChan<Specie>::dpdv
(
    const scalar& rho,
    const scalar& e
) const
{
    scalar V(rho0_/max(rho, 1e-10));
    return
      - Gamma0_*e*sqr(rho)
      - rho*Epsilon1_*A_*pow(V, -Epsilon1_)
      - rho*Epsilon2_*B_*pow(V, -Epsilon2_);
}


template<class Specie>
Foam::scalar Foam::CochranChan<Specie>::dpde
(
    const scalar& rho,
    const scalar& e
) const
{
    return Gamma0_*rho;
}


template<class Specie>
Foam::scalar Foam::CochranChan<Specie>::E
(
    const scalar& rho,
    const scalar& e
) const
{
    scalar V(rho0_/max(rho, 1e-10));
    return
        A_/((Epsilon1_ - 1.0)*rho0_)*pow(V, 1.0 - Epsilon1_)
      - B_/((Epsilon2_ - 1.0)*rho0_)*pow(V, 1.0 - Epsilon2_);
}


template<class Specie>
Foam::scalar Foam::CochranChan<Specie>::Cv
(
    const scalar& rho,
    const scalar& e
) const
{
    return 0;
}


template<class Specie>
Foam::scalar Foam::CochranChan<Specie>::Cp
(
    const scalar& rho,
    const scalar& e
) const
{
    return 0;
}


template<class Specie>
Foam::scalar Foam::CochranChan<Specie>::H
(
    const scalar& rho,
    const scalar& e
) const
{
    scalar V(rho0_/max(rho, 1e-10));
    return
        (
            Epsilon1_*(Epsilon1_ - Gamma0_ - 1.0)*V*pow(V, -Epsilon1_)
          + (Epsilon1_ - 1.0)*Gamma0_*log(rho)
        )*A_/(rho0_*sqr(Epsilon1_ - 1.0))
      + (
            Epsilon2_*(Epsilon2_ - Gamma0_ - 1.0)*V*pow(V, -Epsilon2_)
          + (Epsilon2_ - 1.0)*Gamma0_*log(rho)
        )*B_/(rho0_*sqr(Epsilon2_ - 1.0));
}


template<class Specie>
Foam::scalar Foam::CochranChan<Specie>::CpMCv
(
    const scalar& rho,
    const scalar& e
) const
{
    return this->R();
}


template<class Specie>
Foam::scalar Foam::CochranChan<Specie>::S
(
    const scalar& p,
    const scalar& rho,
    const scalar& T
) const
{
    return -this->R()*log(p/Foam::constant::thermodynamic::Pstd);
}

// ************************************************************************* //
