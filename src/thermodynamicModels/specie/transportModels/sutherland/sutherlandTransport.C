/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "sutherlandTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
Foam::scalar Foam::sutherlandTransport<Thermo>::readCoeff
(
    const word& coeffName,
    const dictionary& dict
)
{
    return dict.subDict("transport").lookup<scalar>(coeffName);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::sutherlandTransport<Thermo>::sutherlandTransport(const dictionary& dict)
:
    Thermo(dict),
    As_(0),
    Ts_(0)
{
    const dictionary& transportDict = dict.subDict("transport");

    if (transportDict.found("As") && transportDict.found("Ts"))
    {
        transportDict.readIfPresent("As", As_);
        transportDict.readIfPresent("Ts", Ts_);
    }
    else if
    (
        transportDict.found("mu0")
     && transportDict.found("T0")
     && transportDict.found("S")
    )
    {
        const scalar mu0 = transportDict.lookup<scalar>("mu0");
        const scalar T0 = transportDict.lookup<scalar>("T0");
        const scalar S = transportDict.lookup<scalar>("S");

        Ts_ = S;
        As_ = mu0/pow(T0, 1.5)*(T0 + S);
    }
    else
    {
        FatalIOErrorInFunction(transportDict)
            << "Missing entries, please provide either" << nl
            << "    As and Ts" << nl
            << "    mu0, T0, S" << endl
            << abort(FatalIOError);
    }
}


template<class Thermo>
Foam::sutherlandTransport<Thermo>::sutherlandTransport
(
    const Thermo& t,
    const dictionary& dict
)
:
    Thermo(t),
    As_(readCoeff("As", dict)),
    Ts_(readCoeff("Ts", dict))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::sutherlandTransport<Thermo>::write(Ostream& os) const
{
    os  << this->name() << endl
        << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    dict.add("As", As_);
    dict.add("Ts", Ts_);

    os  << indent << dict.dictName() << dict
        << decrIndent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const sutherlandTransport<Thermo>& st
)
{
    st.write(os);
    return os;
}


// ************************************************************************* //
