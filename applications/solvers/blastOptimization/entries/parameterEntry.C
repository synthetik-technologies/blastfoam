/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2023-2025
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "parameterEntry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parameterEntry::parameterEntry()
{}



Foam::parameterEntry::parameterEntry(Istream& is)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::parameterEntry::~parameterEntry()
{}


void Foam::parameterEntry::read(Istream& is)
{
    optimizationEntry::read(is);
    is >> value_;
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, parameterEntry& entry)
{
    entry.read(is);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const parameterEntry& entry)
{
    os  << static_cast<const optimizationEntry&>(entry) << incrIndent
        << indent << "value: " << entry.value_ << decrIndent << endl;
    return os;
}


// ************************************************************************* //
