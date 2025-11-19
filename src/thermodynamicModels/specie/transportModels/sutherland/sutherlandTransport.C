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

template<class Thermo>
void Foam::sutherlandTransport<Thermo>::readCoeffs(const dictionary& dict)
{
    const dictionary& transportDict = dict.subDict("transport");

    if
    (
        transportDict.found("Ts")
     || transportDict.found("Smu")
     || transportDict.found("S")
    )
    {
        Smu_ = transportDict.lookupBackwardsCompatible<scalar>
        (
            {"Smu", "Ts", "S"}
        );
    }
    else
    {
        FatalIOErrorInFunction(transportDict)
            << "Could not determine Smu. Please provide, Ts, Smu, or S" << endl
            << abort(FatalIOError);
    }

    if
    (
        transportDict.found("As")
     || transportDict.found("Amu")
    )
    {
        Amu_ = transportDict.lookupBackwardsCompatible<scalar>({"Amu", "As"});
    }
    else if
    (
        transportDict.found("mu0")
     && transportDict.found("T0")
    )
    {
        const scalar T0 = transportDict.lookup<scalar>("T0");
        const scalar mu0 = transportDict.lookup<scalar>("mu0");
        Amu_ = mu0/pow(T0, 1.5)*(T0 + Smu_);
    }
    else
    {
        FatalIOErrorInFunction(transportDict)
            << "Either As/Amu or mu0 and T0 must be provided" << endl
            << abort(FatalIOError);
    }

    // Read kappa information
    if (transportDict.found("Pr"))
    {
        rPr_ = 1.0/transportDict.lookup<scalar>("Pr");
    }
    else if
    (
        transportDict.found("Akappa")
     && transportDict.found("Skappa")
    )
    {
        transportDict.readIfPresent("Akappa", Ak_);
        transportDict.readIfPresent("Skappa", Sk_);
    }
    else if
    (
        transportDict.found("T0")
     && transportDict.found("kappa0")
     && transportDict.found("Skappa")
    )
    {
        const scalar T0 = transportDict.lookup<scalar>("T0");
        const scalar kappa0 = transportDict.lookup<scalar>("kappa0");
        Sk_ = transportDict.lookup<scalar>("Skappa");
        Ak_ = Sk_ > 0 ? kappa0/pow(T0, 1.5)*(T0 + Sk_) : -1;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::sutherlandTransport<Thermo>::sutherlandTransport(const dictionary& dict)
:
    Thermo(dict),
    Amu_(0),
    Smu_(0),
    Ak_(-1),
    Sk_(-1),
    rPr_(-1)
{
    readCoeffs(dict);
}


template<class Thermo>
Foam::sutherlandTransport<Thermo>::sutherlandTransport
(
    const Thermo& t,
    const dictionary& dict
)
:
    Thermo(t),
    Amu_(0),
    Smu_(0),
    Ak_(-1),
    Sk_(-1),
    rPr_(-1)
{
    readCoeffs(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::sutherlandTransport<Thermo>::write(Ostream& os) const
{
    os  << this->name() << endl
        << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    if (Ak_ > 0)
    {
        dict.add("Amu", Amu_);
        dict.add("Smu", Smu_);
        dict.add("Akappa", Ak_);
        dict.add("Skappa", Sk_);
    }
    else
    {
        dict.add("As", Amu_);
        dict.add("Ts", Smu_);
    }

    if (rPr_ > 0)
    {
        dict.add("Pr", 1.0/rPr_);
    }

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
