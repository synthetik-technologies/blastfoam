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
Foam::scalar Foam::Tillotson<Specie>::p
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    if (rho < rho0_)
    {
        return pCav_;
    }

    scalar eTot = e;
    if (rho > rho0_)
    {
        eTot += EcTable_.lookup(rho);
    }

    scalar mu(rho/rho0_ - 1.0);
    return
        p0_
      + omega_*rho*(eTot - e0_)
      + A_*mu
      + B_*sqr(mu)
      + C_*pow3(mu);
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::Pi
(
    const scalar& rho,
    const scalar& e
) const
{
    if (rho < rho0_)
    {
        return 0.0;
    }

    scalar mu(rho/rho0_ - 1.0);
    return
        omega_*rho*e0_
      - p0_
      - A_*mu
      - B_*sqr(mu)
      - C_*pow3(mu);
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::delta
(
    const scalar& p,
    const scalar& rho,
    const scalar& e
) const
{
    if (rho < rhoCav_)
    {
        return 0.0;
    }

    scalar mu(rho/rho0_ - 1.0);
    return
        e0_
      - (
            A_
          + 2.0*B_*mu
          + 3.0*C_*sqr(mu)
        )/rho0_/omega_;
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::speedOfSound
(
    const scalar& p,
    const scalar& rho,
    const scalar& e,
    const scalar& T,
    const scalar& Cv
) const
{
    if (rho < rhoCav_)
    {
        return sqrt(k_*pCav_/max(rho, 1e-10));
    }

    scalar gamma(Gamma(rho, e, T, Cv));
    scalar h((gamma*p + Pi(rho, e))/((gamma - 1.0)*max(rho, 1e-10)));
    scalar cLSqr((h - delta(p, rho, e))*(gamma - 1.0));

    if (rho < rho0_)
    {
        scalar x((rho - rho0_)/(rhoCav_ - rho0_));
        scalar cCavSqr(k_*pCav_/max(rho, 1e-10));

        scalar cSqr
        (
            1.0
           /(
                max(rho, 1e-10)
               *(
                    x/(max(rho, 1e-10)*cCavSqr)
                  + (1.0 - x)/(rho0_*cLSqr)
                )
            )
        );
        return sqrt(max(cSqr, small));
    }
    return sqrt(cLSqr);
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::Gamma
(
    const scalar& rho,
    const scalar& e,
    const scalar& T,
    const scalar& Cv
) const
{
    if (rho < rhoCav_)
    {
        return k_;
    }
    else if (rho < rho0_)
    {
        scalar x((rho - rho0_)/(rhoCav_ - rho0_));

        return
            1.0
           /(
                max(rho, 1e-10)
               *(
                    x/(rhoCav_*k_)
                  + (1.0 - x)/(rho0_*(omega_ + 1.0))
                )
            );
    }

    return omega_ + 1.0;
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::dpdv
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return
        sqr(rho)/pow3(rho0_)
       *(
            pow3(rho0_)*omega_*(e0_ - e)
          - A_*sqr(rho0_)
          - 2.0*B_*rho0_*(rho - rho0_)
          - 3.0*C_*sqr(rho - rho0_)
        );
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::dpde
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return omega_*rho;
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::dpdT
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    NotImplemented;
    return this->R()*rho;
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::E
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::Cv
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    // Check
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::Cp
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    // Check
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::H
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::CpMCv
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return this->R();
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::S
(
    const scalar& p,
    const scalar& rho,
    const scalar& T
) const
{
    return -this->R()*log(p/Foam::constant::thermodynamic::Pstd);
}

// ************************************************************************* //
