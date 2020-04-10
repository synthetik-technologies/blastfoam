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
Foam::scalar Foam::BKW<Specie>::Gamma
(
    const scalar& rho,
    const scalar& T,
    const scalar& cv
) const
{
    return max(dpdT(rho, T)/(cv*max(rho, 1e-10)), small) + 1.0;
}


template<class Specie>
Foam::scalar Foam::BKW<Specie>::p
(
    const scalar& rho,
    const scalar& T
) const
{
    scalar x = kappa_*k_*rho/(this->W()*pow(T + Theta_, a_));
    return this->R()*T*rho*(1.0 + x*exp(beta_*x));
}


template<class Specie>
Foam::scalar Foam::BKW<Specie>::speedOfSound
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
Foam::scalar Foam::BKW<Specie>::dpdv
(
    const scalar& rho,
    const scalar& T
) const
{
    scalar x = kappa_*k_*rho/(this->W()*pow(T + Theta_, a_));
    scalar dxdv = -kappa_*k_*sqr(rho)/(this->W()*pow(T + Theta_, a_));
    return
        this->R()*T*rho
       *(
           exp(beta_*x)*dxdv*(1.0 + x*beta_)
         - (1.0 + x*exp(beta_*x))*rho
        );
}


template<class Specie>
Foam::scalar Foam::BKW<Specie>::dpdT
(
    const scalar& rho,
    const scalar& T
) const
{
    scalar x = kappa_*k_*rho/(this->W()*pow(T + Theta_, a_));
    scalar dxdT = -a_*kappa_*k_*rho/(this->W()*pow(T + Theta_, a_ + 1.0));
    return
        this->R()
       *(
            1.0 + x*exp(beta_*x)
          + T*exp(beta_*x)*dxdT*(1.0 + x*beta_)
        );
}


template<class Specie>
Foam::scalar Foam::BKW<Specie>::E
(
    const scalar& rho,
    const scalar& T
) const
{
    scalar x = kappa_*k_*rho/(this->W()*pow(T + Theta_, a_));
    return this->R()*sqr(T)*a_*x*exp(beta_*x)/(T + Theta_);
}


template<class Specie>
Foam::scalar Foam::BKW<Specie>::Cv
(
    const scalar& rho,
    const scalar& T
) const
{
    scalar x = kappa_*k_*rho/(this->W()*pow(T + Theta_, a_));
    return
        this->R()*T*x*a_*exp(beta_*x)/sqr(T + Theta_)
       *(1.0 + 2.0*Theta_/T - a_*(1.0 + beta_*x));
}


template<class Specie>
Foam::scalar Foam::BKW<Specie>::Cp
(
    const scalar& rho,
    const scalar& T
) const
{
    return Cv(rho, T) - T*sqr(dpdT(rho, T))/dpdv(rho, T) - this->R();
}


template<class Specie>
Foam::scalar Foam::BKW<Specie>::H
(
    const scalar& rho,
    const scalar& T
) const
{
    scalar x = kappa_*k_*rho/(this->W()*pow(T + Theta_, a_));
    return this->R()*T/(T + Theta_)*x*(T*(a_ + 1.0) + Theta_)*exp(beta_*x);
}


template<class Specie>
Foam::scalar Foam::BKW<Specie>::CpMCv
(
    const scalar& rho,
    const scalar& T
) const
{
    scalar x = kappa_*k_*rho/(this->W()*pow(T + Theta_, a_));
    return
        this->R()
       *(
          - T*a_*x/(T + Theta_)*(beta_*x + 1.0)*exp(beta_*x)
          + x*exp(beta_*x) + 1
        );
}


template<class Specie>
Foam::scalar Foam::BKW<Specie>::S
(
    const scalar& p,
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}

// ************************************************************************* //
