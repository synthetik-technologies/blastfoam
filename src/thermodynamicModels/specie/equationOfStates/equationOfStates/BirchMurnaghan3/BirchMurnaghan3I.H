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
Foam::scalar Foam::BirchMurnaghan3<Specie>::p
(
    const scalar& rho,
    const scalar& T
) const
{
    return
        pRef_
      + 1.5*K0_*(pow(rho/rho0_, 7.0/3.0) - pow(rho/rho0_, 5.0/3.0))
       *(1.0 + 0.75*(K0Prime_ - 4.0)*(pow(rho/rho0_, 2.0/3.0) - 1.0));
}


template<class Specie>
Foam::scalar Foam::BirchMurnaghan3<Specie>::Gamma
(
    const scalar& rho,
    const scalar& T,
    const scalar& cv
) const
{
    return Gamma_ + 1.0;
}


template<class Specie>
Foam::scalar Foam::BirchMurnaghan3<Specie>::speedOfSound
(
    const scalar& p,
    const scalar& rho,
    const scalar& T,
    const scalar& cv
) const
{
    scalar cSqr =
        K0_
       *pow(rho/rho0_, 5.0/3.0)
       *(
            6.0*pow(rho/rho0_, 2.0/3.0)*(K0Prime_ - 4.0)
           *(pow(rho/rho0_, 2.0/3.0) - 1.0)
          + (7.0*pow(rho/rho0_, 2.0/3.0) - 5.0)
           *(3.0*(K0Prime_ - 4.0)*(pow(rho/rho0_, 2.0/3.0) - 1.0) + 4.0)
        )/(8.0*max(rho, 1e-10));
    return sqrt(max(cSqr, small));
}


template<class Specie>
Foam::scalar Foam::BirchMurnaghan3<Specie>::dpdv
(
    const scalar& rho,
    const scalar& T
) const
{
    return
      - 0.5*K0_*rho
       *pow(rho/rho0_, 5.0/3.0)
       *(
            6.0*pow(rho/rho0_, 2.0/3.0)*(K0Prime_ - 4.0)
           *(pow(rho/rho0_, 2.0/3.0) - 1.0)
          + (7.0*pow(rho/rho0_, 2.0/3.0) - 5.0)
           *(3.0*(K0Prime_ - 4.0)*(pow(rho/rho0_, 2.0/3.0) - 1.0) + 4.0)
        );
}


template<class Specie>
Foam::scalar Foam::BirchMurnaghan3<Specie>::dpdT
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0;
}


template<class Specie>
Foam::scalar Foam::BirchMurnaghan3<Specie>::E
(
    const scalar& rho,
    const scalar& T
) const
{
    return
      - (
            9.0*K0_*rho/(16.0*pow3(rho0_))
           *(
              - K0Prime_*sqr(rho) + pow(rho, 2.0/3.0)*pow3(rho0_)
               *(
                    3.0*K0Prime_*pow(rho/rho0_, 2.0/3.0)
                  - 3.0*K0Prime_
                  - 14.0*pow(rho/rho0_, 2.0/3.0)
                  + 16.0
                )*pow(rho0_, -5.0/3.0)
              + 4.0*sqr(rho)
            )/16.0
          + pRef_
        )/max(rho, 1e-10);
}


template<class Specie>
Foam::scalar Foam::BirchMurnaghan3<Specie>::Cv
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0;
}


template<class Specie>
Foam::scalar Foam::BirchMurnaghan3<Specie>::H
(
    const scalar& rho,
    const scalar& T
) const
{
    return
        3.0*K0_
       *(
            9.0*K0Prime_*sqr(rho)
          + pow(rho, 2.0/3.0)*pow3(rho0_)
           *(
              - 21.0*K0Prime_*pow(rho/rho0_, 2.0/3.0)
              + 15.0*K0Prime_
              + 98.0*pow(rho/rho0_, 2.0/3.0)
              - 80.0
            )*pow(rho0_, -5.0/3.0)
          - 36.0*sqr(rho)
        )/(16.0*pow3(rho0_));
}


template<class Specie>
Foam::scalar Foam::BirchMurnaghan3<Specie>::Cp
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::BirchMurnaghan3<Specie>::CpMCv
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::BirchMurnaghan3<Specie>::S
(
    const scalar& p,
    const scalar& rho,
    const scalar& T
) const
{
    return 0;
}


// ************************************************************************* //
