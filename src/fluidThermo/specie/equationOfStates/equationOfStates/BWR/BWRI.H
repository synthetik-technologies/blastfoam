/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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
Foam::scalar Foam::BWR<Specie>::Gamma
(
    const scalar& rho,
    const scalar& T,
    const scalar& cv
) const
{
    return dpdT(rho, T)/(cv*max(rho, 1e-10)) + 1.0;
}


template<class Specie>
Foam::scalar Foam::BWR<Specie>::p
(
    const scalar& rho,
    const scalar& T
) const
{
    return
      - this->R()*T*(sqr(rho) + 2.0*pow3(rho)*B0_ + 3.0*pow4(rho)*b_)
      + 2.0*A0_*pow3(rho)
      + 3.0*a_*pow4(rho)
      - 6.0*a_*alpha_*pow(rho, 7)
      + 2.0*C0_*pow3(rho)/sqr(T)
      + c_*pow4(rho)/sqr(T)*exp(-gamma_*sqr(rho))
       *(2.0*sqr(gamma_*sqr(rho)) - 3.0*gamma_*sqr(rho) - 3.0);
}


template<class Specie>
Foam::scalar Foam::BWR<Specie>::speedOfSound
(
    const scalar& p,
    const scalar& rho,
    const scalar& T,
    const scalar& cv
) const
{
    return
        sqrt
        (
            max((sqr(dpdT(rho, T))*T/cv - dpdv(rho, T))/sqr(max(rho, 1e-10)), small)
        );
}


template<class Specie>
Foam::scalar Foam::BWR<Specie>::dpdv
(
    const scalar& rho,
    const scalar& T
) const
{
    return
      - this->R()*T
       *(sqr(rho) + 2.0*pow3(rho)*B0_ + 3.0*pow4(rho)*b_)
      + 2.0*A0_*pow3(rho)
      + 3.0*a_*pow4(rho)
      - 6.0*a_*alpha_*pow(rho, 7)
      + 2.0*C0_*pow3(rho)/sqr(T)
      + c_*pow4(rho)/sqr(T)*exp(-gamma_*sqr(rho))
       *(2.0*sqr(gamma_*sqr(rho)) - 3.0*gamma_*sqr(rho) - 3.0);
}


template<class Specie>
Foam::scalar Foam::BWR<Specie>::dpdT
(
    const scalar& rho,
    const scalar& T
) const
{
    return
        this->R()
       *(
            rho
          + sqr(rho)*B0_
          + pow3(rho)*b_
        )
      + 2.0*C0_*sqr(rho)/pow3(T)
      - 2.0*pow3(rho/T)*c_*(1.0 + gamma_*sqr(rho))*exp(-gamma_*sqr(rho));
}


template<class Specie>
Foam::scalar Foam::BWR<Specie>::E
(
    const scalar& rho,
    const scalar& T
) const
{
    return
      - A0_*rho
      - 3.0*C0_*rho/sqr(T)
      + a_*alpha_*pow5(rho)/5.0
      - a_*sqr(rho)/2.0
      - 3.0*c_/sqr(T)*exp(-gamma_*sqr(rho))*(sqr(rho)/2.0 + 1.0/gamma_);
}


template<class Specie>
Foam::scalar Foam::BWR<Specie>::Cv
(
    const scalar& rho,
    const scalar& T
) const
{
    return
        (6.0*C0_*gamma_*rho + c_*(gamma_*sqr(rho) + 2.0)*exp(-gamma_*sqr(rho)))
       /(pow4(T)*gamma_);
}


template<class Specie>
Foam::scalar Foam::BWR<Specie>::Cp
(
    const scalar& rho,
    const scalar& T
) const
{
    return Cv(rho, T) - T*sqr(dpdT(rho, T))/dpdv(rho, T) - this->R();
}


template<class Specie>
Foam::scalar Foam::BWR<Specie>::H
(
    const scalar& rho,
    const scalar& T
) const
{
    return
        A0_*rho
      + 3.0*C0_*rho/sqr(T)
      - this->R()*T/2.0*(b_*sqr(rho) - 1)
      - 3.0*a_*alpha_*pow5(rho)
      + 1.5*a_*sqr(rho)
      + c_*sqr(rho)*exp(-gamma_*sqr(rho))/sqr(T)
       *(sqr(gamma_)*pow4(rho) - 3.5*(gamma_*sqr(rho) + 1.0));
}


template<class Specie>
Foam::scalar Foam::BWR<Specie>::CpMCv
(
    const scalar& rho,
    const scalar& T
) const
{
    return
        2.0*C0_*rho/pow3(T)
      + this->R()*(B0_*rho + b_*sqr(rho) + 1.0)
      - 2.0*c_*sqr(rho)*(gamma_*sqr(rho) + 1.0)*exp(-gamma_*sqr(rho))/pow3(T);
}


template<class Specie>
Foam::scalar Foam::BWR<Specie>::S
(
    const scalar& p,
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}

// ************************************************************************* //
