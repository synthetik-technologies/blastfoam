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
Foam::scalar Foam::AbelNobel<Specie>::Gamma
(
    const scalar& rho,
    const scalar& T,
    const scalar& cv
) const
{
    return max(dpdT(rho, T)/(cv*max(rho, 1e-10)), small) + 1.0;
}


template<class Specie>
Foam::scalar Foam::AbelNobel<Specie>::p
(
    const scalar& rho,
    const scalar& T
) const
{
    return rho*this->R()*T/(1.0 - b_*rho);
}


template<class Specie>
Foam::scalar Foam::AbelNobel<Specie>::speedOfSound
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
Foam::scalar Foam::AbelNobel<Specie>::dpdv
(
    const scalar& rho,
    const scalar& T
) const
{
    return -this->R()*T*sqr(rho/(1.0 - b_*rho));
}


template<class Specie>
Foam::scalar Foam::AbelNobel<Specie>::dpdT
(
    const scalar& rho,
    const scalar& T
) const
{
    return rho*this->R()/(1.0 - b_*rho);
}


template<class Specie>
Foam::scalar Foam::AbelNobel<Specie>::E
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::AbelNobel<Specie>::Cv
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::AbelNobel<Specie>::Cp
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::AbelNobel<Specie>::H
(
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::AbelNobel<Specie>::CpMCv
(
    const scalar& rho,
    const scalar& T
) const
{
    return this->R();
}


template<class Specie>
Foam::scalar Foam::AbelNobel<Specie>::S
(
    const scalar& p,
    const scalar& rho,
    const scalar& T
) const
{
    return 0.0;
}

// ************************************************************************* //
