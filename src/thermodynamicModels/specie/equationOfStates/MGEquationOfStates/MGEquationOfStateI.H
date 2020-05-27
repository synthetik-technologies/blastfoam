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

template<class EoS>
Foam::scalar Foam::MGEquationOfState<EoS>::p
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return (EoS::Gamma(rho, e) - 1.0)*rho*e - EoS::Pi(rho, e);
}


template<class EoS>
Foam::scalar Foam::MGEquationOfState<EoS>::speedOfSound
(
    const scalar& p,
    const scalar& rho,
    const scalar& e,
    const scalar& T,
    const scalar& Cv
) const
{
    scalar gamma(EoS::Gamma(rho, e));
    scalar h((gamma*p + EoS::Pi(rho, e))/((gamma - 1.0)*max(rho, 1e-10)));
    return sqrt(max((h - EoS::delta(p, rho, e))*(gamma - 1.0), small));
}


template<class EoS>
Foam::scalar Foam::MGEquationOfState<EoS>::Gamma
(
    const scalar& rho,
    const scalar& e,
    const scalar& T,
    const scalar& Cv
) const
{
    return EoS::Gamma(rho, e);
}


template<class EoS>
Foam::scalar Foam::MGEquationOfState<EoS>::dpdv
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return EoS::dpdv(rho, e);
}


template<class EoS>
Foam::scalar Foam::MGEquationOfState<EoS>::dpdT
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    NotImplemented;
    return 0.0;
}


template<class EoS>
Foam::scalar Foam::MGEquationOfState<EoS>::dpde
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return EoS::dpde(rho, e);
}


template<class EoS>
Foam::scalar Foam::MGEquationOfState<EoS>::E
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return EoS::E(rho, e);
}


template<class EoS>
Foam::scalar Foam::MGEquationOfState<EoS>::Cv
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return EoS::Cv(rho, e);
}


template<class EoS>
Foam::scalar Foam::MGEquationOfState<EoS>::Cp
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    // Check
    return EoS::Cp(rho, e);
}


template<class EoS>
Foam::scalar Foam::MGEquationOfState<EoS>::H
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return EoS::H(rho, e);
}


template<class EoS>
Foam::scalar Foam::MGEquationOfState<EoS>::CpMCv
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return EoS::CpMCv(rho, e);
}

// ************************************************************************* //
