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
Foam::scalar Foam::equationOfState<EoS>::p
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return EoS::p(rho, T);
}


template<class EoS>
Foam::scalar Foam::equationOfState<EoS>::speedOfSound
(
    const scalar& p,
    const scalar& rho,
    const scalar& e,
    const scalar& T,
    const scalar& cv
) const
{
    return EoS::speedOfSound(p, rho, T, cv);
}


template<class EoS>
Foam::scalar Foam::equationOfState<EoS>::Gamma
(
    const scalar& rho,
    const scalar& e,
    const scalar& T,
    const scalar& cv
) const
{
    return EoS::Gamma(rho, T, cv);
}


template<class EoS>
Foam::scalar Foam::equationOfState<EoS>::dpdv
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return EoS::dpdv(rho, T);
}


template<class EoS>
Foam::scalar Foam::equationOfState<EoS>::dpdT
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return EoS::dpdT(rho, T);
}


template<class EoS>
Foam::scalar Foam::equationOfState<EoS>::dpde
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    NotImplemented
    return 0;
}


template<class EoS>
Foam::scalar Foam::equationOfState<EoS>::E
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return EoS::E(rho, T);
}


template<class EoS>
Foam::scalar Foam::equationOfState<EoS>::Cv
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return EoS::Cv(rho, T);
}


template<class EoS>
Foam::scalar Foam::equationOfState<EoS>::Cp
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    // Check
    return EoS::Cp(rho, T);
}


template<class EoS>
Foam::scalar Foam::equationOfState<EoS>::H
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return EoS::H(rho, T);
}


template<class EoS>
Foam::scalar Foam::equationOfState<EoS>::CpMCv
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return EoS::CpMCv(rho, T);
}

// ************************************************************************* //
