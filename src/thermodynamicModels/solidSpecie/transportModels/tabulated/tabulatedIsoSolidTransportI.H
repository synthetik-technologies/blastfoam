/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::autoPtr<Foam::tabulatedIsoSolidTransport<Thermo>>
Foam::tabulatedIsoSolidTransport<Thermo>::clone() const
{
    return autoPtr<tabulatedIsoSolidTransport<Thermo>>
    (
        new tabulatedIsoSolidTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::tabulatedIsoSolidTransport<Thermo>>
Foam::tabulatedIsoSolidTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<tabulatedIsoSolidTransport<Thermo>>
    (
        new tabulatedIsoSolidTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::tabulatedIsoSolidTransport<Thermo>::kappa
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return kappaTable_.lookup(T);
}

template<class Thermo>
inline Foam::vector Foam::tabulatedIsoSolidTransport<Thermo>::Kappa
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    const scalar k(kappa(rho, e, T));
    return vector(k, k, k);
}


template<class Thermo>
inline Foam::scalar Foam::tabulatedIsoSolidTransport<Thermo>::mu
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    NotImplemented;
    return scalar(0);
}

// ************************************************************************* //
