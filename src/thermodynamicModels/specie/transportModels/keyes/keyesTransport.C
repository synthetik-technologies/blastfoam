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

#include "keyesTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
Foam::scalar Foam::keyesTransport<Thermo>::readCoeff
(
    const word& coeffName,
    const dictionary& dict
)
{
    return dict.subDict("transport").lookup<scalar>(coeffName);
}

template<class Thermo>
void Foam::keyesTransport<Thermo>::readCoeffs(const dictionary& dict)
{
    const dictionary& transportDict = dict.subDict("transport");

    transportDict.lookup("mu0") >> mu0_;
    transportDict.lookup("B") >> B_;
    transportDict.lookup("C") >> C_;
    rPr_ = 1.0/transportDict.lookup<scalar>("Pr");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::keyesTransport<Thermo>::keyesTransport(const dictionary& dict)
:
    Thermo(dict),
    mu0_(0),
    B_(0),
    C_(1),
    rPr_(-1)
{
    readCoeffs(dict);
}


template<class Thermo>
Foam::keyesTransport<Thermo>::keyesTransport
(
    const Thermo& t,
    const dictionary& dict
)
:
    Thermo(t),
    mu0_(0),
    B_(0),
    C_(1),
    rPr_(-1)
{
    readCoeffs(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::keyesTransport<Thermo>::write(Ostream& os) const
{
    os  << this->name() << endl
        << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    dict.add("mu0", mu0_);
    dict.add("B", B_);
    dict.add("C", C_);
    dict.add("Pr", 1.0/rPr_);

    os  << indent << dict.dictName() << dict
        << decrIndent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const keyesTransport<Thermo>& st
)
{
    st.write(os);
    return os;
}


// ************************************************************************* //
